"""StageSolver -- slice 1: an isentropic (no slip, no losses) impeller.

This is the tracer bullet: the smallest whole machine that proves the
rothalpy-consistent ideal relative total pressure fixes the impeller pressure-ratio
bug (``docs/centrifugal/00-root-cause-analysis.md`` defect A). No diffuser, no slip,
no loss model -- those are added in slices 2-5 through the ``slip=`` / ``losses=``
seams, which this solver already accepts (as ``None``) without changing its shape.

Station numbering follows ``docs/centrifugal/01-design.md`` S3.5: 0 inlet -> 1 impeller
LE -> 2 impeller TE. Slice 1 stops at 2 (no diffuser).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterator, Optional, Sequence

from scipy.optimize import brentq

from .components import Impeller
from .geometry import MeridionalPath
from .losses import (
    NO_LOSSES,
    EvaluatedLosses,
    ImpellerLossState,
    LossSet,
    head_loss_correction,
)
from .state import Air, ThermoState, entropy, rothalpy


class Choked(Exception):
    """Raised when a requested mass flow exceeds what a station can pass at Mm=1.

    Choke is on the MERIDIONAL Mach number (``docs/PHYSICS-RULES.md`` rule 5) -- never the
    absolute Mach, which routinely exceeds 1 at a centrifugal impeller exit without
    choking. Slice 10 needs choke detection anyway; this comes free from bounding the
    continuity solve on ``[eps, V_star]`` (see ``_solve_for_meridional_velocity``).
    """

    def __init__(self, station: str, mdot: float, mdot_choke: float) -> None:
        self.station = station
        self.mdot = mdot
        self.mdot_choke = mdot_choke
        super().__init__(
            f"{station} is choked: requested mdot={mdot:.4f} kg/s exceeds the "
            f"sonic (Mm=1) mass flow {mdot_choke:.4f} kg/s"
        )


@dataclass(frozen=True)
class InletState:
    """Stage inlet boundary condition: absolute total pressure and temperature."""

    P0: float
    T0: float


class StationSeries:
    """The solved stage stations, both iterable and accessible by name."""

    def __init__(self, impeller_le: ThermoState, impeller_te: ThermoState) -> None:
        self.impeller_le = impeller_le
        self.impeller_te = impeller_te

    def __iter__(self) -> Iterator[ThermoState]:
        return iter((self.impeller_le, self.impeller_te))


@dataclass(frozen=True)
class OperatingPoint:
    """The solved operating point: overall performance plus every station state.

    ``losses`` is always an :class:`~turbodesign.centrifugal.losses.EvaluatedLosses`
    (never ``None``) -- empty (``NO_LOSSES``) when ``Stage.losses`` is ``None``, so
    callers never need an ``is None`` check to read ``.internal_total``/``.parasitic_total``.

    ``work_euler`` is the Euler turbomachinery work, ``U2*Vt2`` (no inlet swirl in this
    model). ``work_actual`` is what the shaft actually put into the fluid --
    ``work_euler + losses.parasitic_total``, EXACTLY (docs/PHYSICS-RULES.md rule 6: parasitic ADDS
    work). ``eta_poly`` is the standard polytropic efficiency,
    ``(gamma-1)/gamma * ln(PR) / ln(T02/T01)``, using the ACTUAL (post-loss) PR and
    T-ratio.

    ``impeller_eta_is`` is the isentropic-efficiency split
    ``(work_euler - losses.internal_total) / work_actual``, parasitic in the denominator
    (docs/PHYSICS-RULES.md rule 6). It is **IMPELLER-ONLY, and it is NOT Oh's Eq. (2).**
    This field used to be called ``eta_is`` and to advertise itself as "Oh, Yoon & Chung
    (1997)'s isentropic-efficiency split"; that was wrong twice over.

    1. Oh's Eq. (2) puts SIX internal mechanisms in the numerator, and one of them is the
       VANELESS-DIFFUSER loss -- a STATIONARY component. Our ``internal_total``
       (``Stage.solve``) sums only the SEVEN IMPELLER internal models. Same algebraic
       shape, different loss set: it is not Oh's number and must not be reported as such.
    2. The old docstring's parenthetical "``dh_static = 0`` -- no stationary-component
       losses exist yet; ``diffuser`` is still ``None``-only" went STALE the moment
       ``Stage.components`` was populated. Stationary losses DO exist now; this quantity
       simply does not include them, BY CONSTRUCTION.

    For a number comparable to a measured stage efficiency, use ``stage_eta_poly_realgas``
    (post-diffuser, post-EGV, real-gas exponent) -- that is what the paper reports.
    """

    PR: float  # IMPELLER total-to-total pressure ratio (P02/P01) -- NOT the stage value
    psi: float  # work factor, Cp*dT0 / U2**2
    stations: StationSeries
    losses: EvaluatedLosses = NO_LOSSES
    work_euler: float = 0.0
    work_actual: float = 0.0
    eta_poly: float = 1.0  # IMPELLER polytropic efficiency -- NOT the stage value
    # IMPELLER isentropic efficiency -- NOT the stage value, and NOT Oh's Eq. (2) (which
    # also carries the vaneless-diffuser loss in its numerator). See the class docstring.
    impeller_eta_is: float = 1.0

    # The stationary flow path, when Stage.components is configured. THESE are the numbers
    # comparable to NASA's published STAGE values (PR_tt 4.6847, eta_poly 0.8553);
    # ``PR``/``eta_poly`` above are impeller-only and are comparable to NOTHING. Keeping
    # both, separately named, so no one can accidentally compare the wrong pair -- that
    # confusion is precisely what let psi be the sole constrained output.
    stage_states: tuple = ()
    stage_PR: Optional[float] = None  # P0_exit / P01
    stage_eta_poly: Optional[float] = None  # constant-cp exponent -- INFLATED ~0.8 pts
    # Real-gas eta: polytropic exponent with cp mass-averaged over the actual compression.
    # THIS is the value comparable to NASA's measurement. See Air.cp_at.
    stage_eta_poly_realgas: Optional[float] = None

    @property
    def stage_exit(self):
        """The final stationary state, or None if no components are configured."""
        return self.stage_states[-1] if self.stage_states else None


def _sonic_root(C0: float, C1: float, C2: float, fluid: Air) -> float:
    """Positive root of ``Cm**2 = gamma*R*(C0 + C1*Cm - C2*Cm**2)`` -- the meridional
    velocity where the MERIDIONAL Mach number equals 1, for a local static temperature
    ``T(Cm) = C0 + C1*Cm - C2*Cm**2``.

    This generalises :func:`_sonic_meridional_velocity` to a T(Cm) that has a LINEAR
    term as well as the quadratic one. That linear term is exactly zero in slice 1 (no
    slip: ``Vt2 = U2 + Cm2*tan(beta2b)``, and the linear-in-Cm2 pieces cancel -- see
    that function's docstring), which is why slice 1 never needed this. Slip breaks
    the cancellation: ``Vt2 = sigma*U2 + Cm2*tan(beta2b)`` with ``sigma < 1`` decouples
    the constant offset from ``U2``, so ``T2(Cm2)`` picks up a genuine linear term (see
    the derivation in ``Stage.solve``'s impeller-exit block). Substituting
    ``Cm**2 = gamma*R*T(Cm)`` gives a quadratic in ``Cm``, solved directly below --
    still closed-form, no iteration.

    Reduces exactly to ``_sonic_meridional_velocity(T0_eff, k, fluid)`` when ``C1 = 0``
    (verified: with ``a = 1 + gamma*R*C2``, ``b = 0``, the positive root is
    ``sqrt(gamma*R*C0/a)``, which is ``_sonic_meridional_velocity``'s formula with
    ``T0_eff = C0`` and ``(gamma-1)*k/2 = gamma*R*C2``).
    """
    gR = fluid.gamma * fluid.R
    a = 1.0 + gR * C2
    b = -gR * C1
    c = -gR * C0
    disc = b * b - 4.0 * a * c
    return (-b + math.sqrt(disc)) / (2.0 * a)


def _sonic_meridional_velocity(T0_eff: float, k: float, fluid: Air) -> float:
    """The meridional velocity ``V_star`` where the MERIDIONAL Mach number equals 1.

    Valid at any station whose local static temperature can be written
    ``T(Vm) = T0_eff - k*Vm**2/(2*cp)``: ``k = 1`` at an inlet with no swirl (the
    absolute velocity is purely meridional, so this is exactly ``T1 = T01 -
    Cx1**2/(2*cp)``); ``k = 1 + tan(beta2b)**2`` at the impeller exit, because
    ``Vt2 = U2 + Cm2*tan(beta2b)`` couples the swirl -- and hence T02 and the
    ``Vt2**2`` term -- back to Cm2, and the linear-in-Cm2 terms cancel exactly,
    leaving ``T2(Cm2) = T01 + U2**2/(2*cp) - (1+tan(beta2b)**2)*Cm2**2/(2*cp)``, the
    same shape with ``T0_eff = T01 + U2**2/(2*cp)``.

    Substituting into ``Vm**2 = gamma*R*T(Vm)`` and solving for ``Vm`` (using
    ``gamma*R = cp*(gamma-1)``) gives the closed form below -- no iteration needed.
    """
    R = fluid.R
    gamma = fluid.gamma
    return math.sqrt(gamma * R * T0_eff / (1.0 + (gamma - 1.0) * k / 2.0))


def _solve_for_meridional_velocity(
    residual, V_star: float, station: str, mdot: float, eps: float = 1e-6
):
    """Solve ``residual(Vm) = 0`` for Vm on ``(0, V_star]`` -- the subsonic branch only.

    ``V_star`` is the sonic (Mm=1) meridional velocity (see
    ``_sonic_meridional_velocity``). The mass-flux residual is monotone increasing on
    this interval, so it brackets a root IFF ``mdot <= mdot_choke`` -- no outward grid
    walk, no risk of stepping over both the subsonic and supersonic roots into
    negative-temperature territory (docs/PHYSICS-RULES.md rule 5: bound on the MERIDIONAL Mach,
    never on ``|V|`` -- at the impeller exit M2 can be ~1.08 while Mm2 is ~0.22, so
    bounding on |V| would reject perfectly valid, unchoked operating points).

    If ``residual(V_star) < 0``, the station cannot pass the requested mass flow even
    at Mm=1: raise :class:`Choked` by name instead of walking on into ``T < 0``, where
    ``(negative)**n`` for non-integer ``n`` returns a complex number.
    """
    f_lo = residual(eps)
    f_hi = residual(V_star)
    if f_hi < 0.0:
        raise Choked(station, mdot, mdot_choke=f_hi + mdot)
    if f_lo >= 0.0:
        # mdot ~ 0 (or negative): no physical operating point on this branch either.
        raise Choked(station, mdot, mdot_choke=f_lo + mdot)
    return brentq(residual, eps, V_star, xtol=1e-10)


class NotConverged(Exception):
    """Raised when a coupled residual solve fails to converge -- never returned
    silently (an unconverged answer is worse than a loud failure).
    """


def _solve_coupled_meridional_velocity(
    residual,
    V_star: float,
    station: str,
    mdot: float,
    eps: float = 1e-6,
    max_iter: int = 200,
    residual_tol: float = 1e-8,
):
    """Like :func:`_solve_for_meridional_velocity`, but for a ``residual`` that folds
    a Cm-DEPENDENT correction (e.g. an internal loss's effect on rho2) into the mass-
    flux balance -- i.e. a genuinely coupled continuity solve, not a pure lossless
    mass-flux balance (see the impeller-exit block of ``Stage.solve``).

    ``brentq``'s ``xtol`` bounds the bracket on Cm, not the physical mass-flux
    residual it implies once density itself is a function of the trial Cm. This
    wrapper additionally verifies Brent's own convergence flag AND re-evaluates the
    residual at the returned root, checking it is below ``residual_tol`` (relative to
    mdot) within ``max_iter`` iterations -- RAISING :class:`NotConverged` rather than
    ever returning an answer that only *looks* converged from Cm's bracket width.

    Returns ``(Cm, iterations)``.
    """
    f_lo = residual(eps)
    f_hi = residual(V_star)
    if f_hi < 0.0:
        raise Choked(station, mdot, mdot_choke=f_hi + mdot)
    if f_lo >= 0.0:
        raise Choked(station, mdot, mdot_choke=f_lo + mdot)
    Cm, info = brentq(residual, eps, V_star, xtol=1e-10, maxiter=max_iter, full_output=True)
    if not info.converged:
        raise NotConverged(
            f"{station}: coupled continuity solve did not converge within {max_iter} iterations"
        )
    rel_resid = abs(residual(Cm)) / mdot
    if rel_resid > residual_tol:
        raise NotConverged(
            f"{station}: coupled continuity residual {rel_resid:.3e} (relative to "
            f"mdot) exceeds tolerance {residual_tol:.1e} after {info.iterations} "
            "Brent iterations"
        )
    return Cm, info.iterations


@dataclass(frozen=True)
class Stage:
    """A centrifugal compressor stage: geometry + components + fluid.

    ``slip`` accepts any ``turbodesign.centrifugal.slip.SlipModel`` (``None`` means
    perfect blade guidance, slice 1's assumption). ``losses`` accepts a
    ``turbodesign.centrifugal.losses.LossSet`` (``None`` means a lossless impeller,
    slices 1/2's assumption) -- see :mod:`turbodesign.centrifugal.losses` for the
    internal/parasitic split and how each is applied.
    ``diffuser`` is likewise accepted but must be ``None`` here: silently ignoring a
    configured diffuser would return the impeller-only PR with no error and no
    warning, exactly the ``LossType.Enthalpy`` failure mode (a configured component
    that contributes nothing -- ``docs/PHYSICS-RULES.md``'s "single most important thing to know"
    table). ``op.PR`` is the impeller-only pressure ratio until slice 3 fills this in.
    """

    path: MeridionalPath
    impeller: Impeller
    diffuser: Optional[object] = None
    slip: Optional[object] = None
    losses: Optional[LossSet] = None
    fluid: Air = Air()
    # The stationary flow path downstream of the impeller, solved in order: vaneless
    # space -> vaned diffuser -> EGV (turbodesign.centrifugal.diffusion). Empty means the
    # model stops at the impeller TE, in which case op.PR and op.eta_poly are
    # IMPELLER-ONLY and are NOT comparable to NASA's stage values -- which is precisely
    # what left psi as the only constrained output and made the fit unidentifiable
    # (docs/centrifugal/07-diffusion-system-design.md).
    components: Sequence[object] = ()

    def solve(self, mdot: float, rpm: float, inlet: InletState) -> OperatingPoint:
        if self.diffuser is not None:
            raise NotImplementedError(
                "the legacy `diffuser=` slot is not implemented -- pass the stationary "
                "flow path as `components=[VanelessSpace(...), VanedDiffuser(...), "
                "ExitGuideVanes(...)]` instead"
            )

        fluid = self.fluid
        cp, R = fluid.cp, fluid.R
        # Isentropic exponent as cp/R, not gamma/(gamma-1). These are equal exactly:
        # state.Air derives gamma from (cp, R) precisely so cp/R == gamma/(gamma-1) to
        # machine precision (see Air's docstring). Written as cp/R here only because
        # that is the form the entropy relation below needs directly -- not because it
        # differs from gamma/(gamma-1). Do not reintroduce an independently-supplied
        # gamma; that is exactly what breaks the identity (Air's docstring works the
        # cp=1005/gamma=1.4/R=287 example that manufactured ds=+0.185 J/(kg K) across a
        # lossless impeller).
        n = cp / R
        omega = rpm * math.pi / 30.0
        T01, P01 = inlet.T0, inlet.P0

        # ---- impeller inlet (the inducer eye; see geometry.MeridionalPath.inducer_eye)
        le_geom = (
            self.path.station_at_x(self.impeller.x_le)
            if self.impeller.x_le is not None
            else self.path.inducer_eye()
        )
        r1h, r1s = le_geom.r_hub, le_geom.r_shroud
        A1 = le_geom.area  # blockage=0.0 by default: no boundary-layer assumption at inlet
        # Streamline-consistent: the whole relative-frame closure is evaluated on ONE
        # streamline, the RMS radius (area-representative for an annulus).
        r1_rms = math.sqrt(0.5 * (r1h**2 + r1s**2))
        U1 = omega * r1_rms

        def le_residual(Cx: float) -> float:
            T1 = T01 - Cx**2 / (2.0 * cp)
            if T1 <= 0.0:
                return -mdot  # guard: never let a negative T reach T1**n (complex result)
            P1 = P01 * (T1 / T01) ** n
            rho1 = P1 / (R * T1)
            return rho1 * Cx * A1 - mdot

        V_star1 = _sonic_meridional_velocity(T0_eff=T01, k=1.0, fluid=fluid)
        Cx1 = _solve_for_meridional_velocity(le_residual, V_star1, "impeller inlet", mdot)
        T1 = T01 - Cx1**2 / (2.0 * cp)
        P1 = P01 * (T1 / T01) ** n
        rho1 = P1 / (R * T1)
        assert Cx1 < fluid.speed_of_sound(T1), (
            "impeller inlet root is not on the subsonic branch (Mm >= 1) -- the "
            "[eps, V_star] bound should make this unreachable; something upstream "
            "moved the bracket"
        )
        Vm1 = Cx1
        Vx1 = Vm1 * math.cos(le_geom.phi)
        Vr1 = Vm1 * math.sin(le_geom.phi)
        W1 = math.hypot(Vm1, U1)  # no inlet prewhirl: Vt1 = 0
        T0R1 = T1 + W1**2 / (2.0 * cp)
        P0R1 = P1 * (T0R1 / T1) ** n
        I1 = rothalpy(T1, W1, U1, fluid)
        s1 = entropy(T1, P1, T01, P01, fluid)

        le_state = ThermoState(
            P0=P01,
            T0=T01,
            P=P1,
            T=T1,
            P0R=P0R1,
            P0R_is=P0R1,
            T0R=T0R1,
            rothalpy=I1,
            U=U1,
            Vm=Vm1,
            Vx=Vx1,
            Vr=Vr1,
            Vt=0.0,
            W=W1,
            s=s1,
        )

        # ---- impeller exit
        te_geom = self.path.station_at_radius(self.impeller.r_te).with_blockage(
            self.impeller.blockage
        )
        A2 = te_geom.area
        U2 = omega * self.impeller.r_te
        tan_b2b = math.tan(math.radians(self.impeller.backsweep_deg))
        Z_eff = self.impeller.Z_eff

        def _vt2(Cm2: float) -> float:
            """Vt2 = sigma*U2 + Cm2*tan(beta2b) (docs/PHYSICS-RULES.md rule 7): a velocity
            deficit, never an additive angle. ``sigma = 1`` (perfect blade guidance)
            when no slip model is configured -- slice 1's assumption, preserved
            exactly as the ``self.slip is None`` branch.

            ``dbeta_dm``/``r2`` (docs/centrifugal/17-tdd-plan.md slice S6) are passed
            ONLY through the ``SlipModel`` protocol's ``**kw`` seam: every model
            except ``QiuSlip`` ignores them via its own ``**kw``, so this is
            output-INVARIANT for ``WiesnerSlip`` -- the default in every fixture
            (``tests/unit/test_qiu_slip.py::test_shipped_outputs_are_bit_identical``).
            """
            if self.slip is None:
                return U2 + Cm2 * tan_b2b
            sigma = self.slip.sigma(
                beta2b_deg=self.impeller.backsweep_deg,
                Z=Z_eff,
                Cm2_over_U2=Cm2 / U2,
                omega=omega,
                dbeta_dm=self.impeller.dbeta_dm_te,
                r2=self.impeller.r_te,
            )
            return sigma * U2 + Cm2 * tan_b2b

        def _internal_losses_at(
            Cm2: float, Vt2: float, W2: float, T2: float, rho2_is: float
        ) -> dict:
            """Every configured INTERNAL loss's delta_h at a trial (Cm2, Vt2, W2, T2)
            velocity triangle, keyed by class name.

            COUPLED into the continuity solve (``te_residual`` below calls this on
            every trial Cm2): an internal loss destroys P02, which lowers rho2, which
            -- at fixed mdot -- RAISES the Cm2 continuity settles on. Because
            ``Vt2 = sigma*U2 + Cm2*tan(beta2b)`` with ``beta2b < 0`` (backswept), a
            higher Cm2 then LOWERS Vt2, and hence the Euler work and T02 itself. That
            is real physics, not noise: at the HECC design point it is a ~2% P02 drop,
            an ~86.9->88.8 m/s Cm2 rise, and a ~0.7 K T02 fall. Suppressing the
            coupling (evaluating this only after a lossless continuity solve, as a
            prior version of this module did) makes rho2/Cm2 self-inconsistent with
              the loss actually applied, and the error grows off-design.

            What internal loss still does NOT do is add a work term: the Euler-work
            formula ``w = U2*Vt2`` never sees ``internal_total`` directly -- only
            through Vt2's dependence on the converged Cm2. That is the exact,
            unconditional invariant this slice's tests check
            (``test_internal_loss_adds_no_work_term``): ``work_actual == work_euler``
            bit-for-bit whenever no parasitic loss is configured.

            ``rho2_is`` is the TRIAL ISENTROPIC (loss-free) exit density at this Cm2 --
            the caller already holds it as an intermediate (from ``P0R2_is``) before
            applying any internal loss, so passing it through costs nothing and lets
            an internal loss that itself needs a rho2/rho1 compressibility ratio
            (Jansen clearance, slice 5) avoid depending circularly on its own output.
            See ``turbodesign.centrifugal.losses.ImpellerLossState``'s docstring.
            """
            if self.losses is None:
                return {}
            state = ImpellerLossState(
                fluid=fluid,
                mdot=mdot,
                omega=omega,
                Z=Z_eff,
                # MAIN blades only at the LE -- the splitters start downstream and cannot
                # block the inducer eye (ImpellerLossState.Z_le).
                Z_le=float(self.impeller.n_blades),
                backsweep_deg=self.impeller.backsweep_deg,
                r1h=r1h,
                r1s=r1s,
                U1=U1,
                Vm1=Vm1,
                T1=T1,
                rho1=rho1,
                r2=self.impeller.r_te,
                b2=te_geom.b,
                tip_clearance=self.impeller.tip_clearance,
                beta1b_deg=self.impeller.inducer_blade_angle_deg,
                t1=self.impeller.le_blade_thickness,
                A_th=self.impeller.throat_area,
                # slice S1 (docs/centrifugal/17-tdd-plan.md): the two DIFFERENT length
                # quantities the camberline-bug fix separated. L_blade (through-blade)
                # feeds skin friction + mixing; L_meridional feeds leakage only. Both
                # None -> ImpellerLossState's own fallbacks (the corrected closure).
                L_blade=self.impeller.l_blade_m,
                L_meridional=self.impeller.l_main_m,
                # slice S2 (docs/centrifugal/17-tdd-plan.md): the diffusion factor's
                # OWN blade count/K_BL, sourced from Galvas (Eq. B59 + FORTRAN
                # CONST1/SPLT branch) -- an INTEGER exit count (both blade rows reach
                # the TE), NOT Z_eff. n_splitters == 0 -> Z_exit == n_blades == Z_eff
                # and has_splitters == False, so an unsplittered impeller (Eckardt
                # O/A) is bit-identical (ImpellerLossState.Z_exit's docstring).
                Z_exit=float(self.impeller.n_blades + self.impeller.n_splitters),
                has_splitters=self.impeller.n_splitters > 0,
                # slice S4 (docs/centrifugal/17-tdd-plan.md): the SAME blockage the
                # velocity triangle already applied to te_geom.area above
                # (`te_geom = self.path.station_at_radius(...).with_blockage(
                # self.impeller.blockage)`). Without this, ImpellerMixingAungier's own
                # A2 = 2*pi*r2*b2 is UNBLOCKED while Cm2 was solved against a BLOCKED
                # area -- the exact silent decoupling this slice exists to close. Site
                # ONE of TWO -- the parasitic-block construction below is the other;
                # missing either is the bug class this project keeps finding (Z_le,
                # then Z_exit/has_splitters, now this).
                blockage2=self.impeller.blockage,
                U2=U2,
                Cm2=Cm2,
                Vt2=Vt2,
                W2=W2,
                T2=T2,
                rho2=rho2_is,
            )
            internal = {
                type(m).__name__: m.delta_h(state)
                for m in self.losses.models
                if m.kind == "internal"
            }
            # Aungier (1995) Eq. (33): mu_2T = I_B - f_c*SUM(dq). f_c scales the INTERNAL
            # sum ONLY -- the parasitic terms are separate work inputs in his Eq. (1), and
            # scaling them would ADD WORK no source puts there (rule 6). ADOPTED E-R11
            # (C3); see losses.head_loss_correction for the primary text and the
            # hybrid-set assumption this rests on.
            #
            # *** APPLIED HERE, NOT IN LossSet.evaluate. *** This function -- NOT
            # evaluate() -- is what feeds ds_internal and therefore the pressure solve;
            # evaluate() is called once at the end, for REPORTING. f_c placed only in
            # evaluate() is a TERM THAT CANNOT FIRE: it would change the reported loss
            # breakdown and NOTHING ELSE, silently leaving PR untouched. That is exactly
            # how this codebase's dead-loss defects have always looked.
            if getattr(self.losses, "apply_head_loss_correction", False) and internal:
                fc = head_loss_correction(state)
                internal = {k: fc * v for k, v in internal.items()}
            return internal

        def te_residual(Cm2: float) -> float:
            """Mass-flux residual at trial Cm2, COUPLED to internal loss.

            The Euler-work-only velocity triangle (w, T02, W2, T2) never sees the
            loss -- internal loss adds no work term, ever. But the DENSITY used for
            continuity does: internal loss's ds = dh/T2 (Denton 1993 S2) knocks P0R2,
            and hence P2 and rho2, down from the loss-free (``_is``) value at THIS
            trial Cm2 -- so raising Cm2 to compensate for the lower density is exactly
            what the root-finder does by construction. Parasitic loss never appears
            here (docs/PHYSICS-RULES.md rule 6: it adds work, not pressure -- it must not touch
            continuity/density; see the module docstring of
            ``turbodesign.centrifugal.losses``).
            """
            Vt2 = _vt2(Cm2)
            w = U2 * Vt2  # Euler work, no inlet swirl
            T02 = T01 + w / cp
            W2 = math.hypot(Cm2, Vt2 - U2)
            T2 = T02 - (Cm2**2 + Vt2**2) / (2.0 * cp)
            if T2 <= 0.0:
                return -mdot  # guard: never let a negative T reach T2**n (complex result)
            T0R2 = T2 + W2**2 / (2.0 * cp)
            # THE fix: never P0R2_is = P0R1 (valid only if U1 == U2).
            P0R2_is = P0R1 * (T0R2 / T0R1) ** n
            # Trial ISENTROPIC (loss-free) exit density at this Cm2 -- an intermediate
            # the internal-loss/continuity solve needs anyway, exposed to
            # _internal_losses_at so a loss with its own rho2/rho1 dependence (Jansen
            # clearance) does not become circular with the internal-loss total it is
            # part of (see that function's docstring).
            P2_is = P0R2_is * (T2 / T0R2) ** n
            rho2_is = P2_is / (R * T2)
            internal = _internal_losses_at(Cm2, Vt2, W2, T2, rho2_is)
            ds_internal = sum(internal.values()) / T2
            P0R2 = P0R2_is * math.exp(-ds_internal / R)  # Denton 1993 S2
            P2 = P0R2 * (T2 / T0R2) ** n
            rho2 = P2 / (R * T2)
            return rho2 * Cm2 * A2 - mdot

        if self.slip is None:
            # T2(Cm2) = T0_eff - (1+tan(beta2b)**2)*Cm2**2/(2cp): see
            # _sonic_meridional_velocity's docstring for the derivation of T0_eff and k.
            T0_eff2 = T01 + U2**2 / (2.0 * cp)
            V_star2 = _sonic_meridional_velocity(T0_eff2, 1.0 + tan_b2b**2, fluid)
        else:
            # Slip decouples the tangential-velocity offset from U2 (sigma < 1), so the
            # linear-in-Cm2 terms that cancelled exactly above no longer do: T2(Cm2)
            # becomes a genuine quadratic C0 + C1*Cm2 - C2*Cm2**2 (see _sonic_root's
            # docstring for the derivation). sigma is evaluated once here at a
            # reference operating point purely to get a search BRACKET -- exact for
            # the geometry-only models (Wiesner/Stanitz/Busemann, whose sigma does not
            # depend on Cm2_over_U2 at all) and an adequate approximation for Qiu's
            # flow-dependent sigma, because te_residual/_vt2 re-evaluate sigma at the
            # true trial Cm2 on every call -- the bracket only needs to contain the
            # root, not equal it.
            sigma_ref = self.slip.sigma(
                beta2b_deg=self.impeller.backsweep_deg,
                Z=Z_eff,
                Cm2_over_U2=None,
                omega=omega,
                dbeta_dm=self.impeller.dbeta_dm_te,
                r2=self.impeller.r_te,
            )
            a = sigma_ref * U2
            C0 = T01 + U2 * a / cp - a**2 / (2.0 * cp)
            C1 = tan_b2b * (U2 - a) / cp
            C2 = (1.0 + tan_b2b**2) / (2.0 * cp)
            V_star2 = _sonic_root(C0, C1, C2, fluid)
        # Coupled continuity/internal-loss solve (Tier B: residual < 1e-8 relative to
        # mdot, <= 200 Brent iterations) -- RAISES rather than ever returning an
        # unconverged answer. ``_n_te_iterations`` is not consumed further; it exists
        # so this is inspectable (e.g. in a debugger/REPL) without changing the
        # return shape of Stage.solve.
        Cm2, _n_te_iterations = _solve_coupled_meridional_velocity(
            te_residual, V_star2, "impeller exit", mdot
        )
        Vt2 = _vt2(Cm2)
        w = U2 * Vt2  # Euler work, no inlet swirl -- work_euler below
        T02_aero = T01 + w / cp  # Euler-work-only T02; the parasitic bump is added below
        W2 = math.hypot(Cm2, Vt2 - U2)
        T2 = T02_aero - (Cm2**2 + Vt2**2) / (2.0 * cp)
        assert Cm2 < fluid.speed_of_sound(T2), (
            "impeller exit root is not on the subsonic branch (Mm >= 1) -- the "
            "[eps, V_star] bound should make this unreachable; something upstream "
            "moved the bracket"
        )
        T0R2 = T2 + W2**2 / (2.0 * cp)
        P0R2_is = P0R1 * (T0R2 / T0R1) ** n
        P2_is = P0R2_is * (T2 / T0R2) ** n
        rho2_is = P2_is / (R * T2)

        # Recompute the internal-loss/pressure state at the CONVERGED Cm2 -- bit-for-
        # bit the same values ``te_residual`` used on its last (converged) call, just
        # named for the rest of this function (rho2 now feeds the parasitic block and
        # the returned station state).
        internal = _internal_losses_at(Cm2, Vt2, W2, T2, rho2_is)
        internal_total = sum(internal.values())
        ds_internal = internal_total / T2  # same conversion as in te_residual, at Cm2's root
        P0R2 = P0R2_is * math.exp(-ds_internal / R)
        P2 = P0R2 * (T2 / T0R2) ** n
        rho2 = P2 / (R * T2)
        # Internal loss only: P02 uses T02_aero, never the parasitic-bumped T02 below --
        # parasitic work adds NO pressure (docs/PHYSICS-RULES.md rule 6), so P02 must be bit-for-bit
        # independent of which parasitic models are configured.
        P02 = P2 * (T02_aero / T2) ** n

        parasitic: dict = {}
        if self.losses is not None:
            state2 = ImpellerLossState(
                fluid=fluid,
                mdot=mdot,
                omega=omega,
                Z=Z_eff,
                # MAIN blades only at the LE -- the splitters start downstream and cannot
                # block the inducer eye (ImpellerLossState.Z_le).
                Z_le=float(self.impeller.n_blades),
                backsweep_deg=self.impeller.backsweep_deg,
                r1h=r1h,
                r1s=r1s,
                U1=U1,
                Vm1=Vm1,
                T1=T1,
                rho1=rho1,
                r2=self.impeller.r_te,
                b2=te_geom.b,
                tip_clearance=self.impeller.tip_clearance,
                beta1b_deg=self.impeller.inducer_blade_angle_deg,
                t1=self.impeller.le_blade_thickness,
                A_th=self.impeller.throat_area,
                # slice S1 (docs/centrifugal/17-tdd-plan.md): the two DIFFERENT length
                # quantities the camberline-bug fix separated. L_blade (through-blade)
                # feeds skin friction + mixing; L_meridional feeds leakage only. Both
                # None -> ImpellerLossState's own fallbacks (the corrected closure).
                L_blade=self.impeller.l_blade_m,
                L_meridional=self.impeller.l_main_m,
                # slice S2 (docs/centrifugal/17-tdd-plan.md): SAME reasoning as the
                # internal-loss construction site above -- this is the SECOND of the
                # two ImpellerLossState sites this slice must populate. Missing one
                # would silently split blade-loading (uses the first site's state via
                # the coupled trial) and recirculation (uses THIS state) onto two
                # different D_f's, which is exactly the bug class this project keeps
                # finding (Z_le landed on only one site once too).
                Z_exit=float(self.impeller.n_blades + self.impeller.n_splitters),
                has_splitters=self.impeller.n_splitters > 0,
                # slice S4 (docs/centrifugal/17-tdd-plan.md): SITE TWO of TWO -- see the
                # internal-loss construction site's comment above for why both must be
                # populated (ImpellerMixingAungier is internal-only today and reads
                # only the first site's state, but the SAME two-site discipline this
                # project already enforces for Z_exit/has_splitters applies here: a
                # field silently missing from one site is exactly how Z_le landed on
                # only one once).
                blockage2=self.impeller.blockage,
                U2=U2,
                Cm2=Cm2,
                Vt2=Vt2,
                W2=W2,
                T2=T2,
                rho2=rho2,
            )
            for m in self.losses.models:
                if m.kind == "parasitic":
                    parasitic[type(m).__name__] = m.delta_h(state2)
                elif m.kind != "internal":
                    raise ValueError(
                        f"{type(m).__name__}.kind={m.kind!r} is neither 'internal' nor "
                        "'parasitic' -- LossModel.kind has no default because guessing "
                        "wrong inverts the physics (docs/PHYSICS-RULES.md rule 6)."
                    )
        parasitic_total = sum(parasitic.values())
        losses_result = (
            EvaluatedLosses(internal=internal, parasitic=parasitic)
            if self.losses is not None
            else NO_LOSSES
        )

        # Parasitic loss ADDS work with NO pressure benefit (docs/PHYSICS-RULES.md rule 6): it
        # raises the reported T02 only. rho2/P2/Cm2/W2/T2 above -- and hence rothalpy
        # and entropy below -- stay on the aerodynamic (parasitic-free) state: the
        # disc-friction/windage heat this bookkeeps is generated on the impeller BACK
        # FACE, outside the primary flow passage this meanline solves for.
        T02 = T02_aero + parasitic_total / cp
        Vx2 = Cm2 * math.cos(te_geom.phi)
        Vr2 = Cm2 * math.sin(te_geom.phi)
        I2 = rothalpy(T2, W2, U2, fluid)
        s2 = entropy(T2, P2, T01, P01, fluid)

        te_state = ThermoState(
            P0=P02,
            T0=T02,
            P=P2,
            T=T2,
            P0R=P0R2,
            P0R_is=P0R2_is,
            T0R=T0R2,
            rothalpy=I2,
            U=U2,
            Vm=Cm2,
            Vx=Vx2,
            Vr=Vr2,
            Vt=Vt2,
            W=W2,
            s=s2,
        )

        PR = P02 / P01
        psi = cp * (T02 - T01) / U2**2

        work_euler = w
        work_actual = work_euler + parasitic_total  # exact, Tier A (rel<=1e-9)
        gamma = fluid.gamma
        eta_poly = ((gamma - 1.0) / gamma) * math.log(PR) / math.log(T02 / T01)
        # IMPELLER isentropic efficiency: parasitic in the denominator, internal in the
        # numerator (docs/PHYSICS-RULES.md rule 6). internal_total sums the SEVEN IMPELLER
        # internal models ONLY -- the stationary components' losses are NOT in it, so this
        # is NOT Oh 1997's Eq. (2) (his numerator carries the vaneless-diffuser loss too).
        # Compare a measured stage efficiency against stage_eta_poly_realgas, not this.
        impeller_eta_is = (work_euler - internal_total) / work_actual

        # ---------------------------------------------------- the stationary flow path
        # Solve the components in order. T0 is conserved across all of them (U = 0 -> no
        # work), so the STAGE work factor is the impeller's -- but the stage PR and eta
        # are not, because each component destroys P0. These are the numbers that are
        # actually comparable to NASA's published stage values.
        stage_states: list[object] = []
        stage_exit = None
        if self.components:
            from .diffusion import state_from_totals

            cur = state_from_totals(
                P0=te_state.P0,
                T0=te_state.T0,
                Vm=te_state.Vm,
                Vt=te_state.Vt,
                r=self.impeller.r_te,
                b=te_geom.b,
                fluid=fluid,
                s=te_state.s,
            )
            for comp in self.components:
                cur = comp.solve(cur, mdot, fluid)
                stage_states.append(cur)
            stage_exit = cur

        stage_PR = stage_exit.P0 / P01 if stage_exit is not None else None
        stage_eta_poly = None
        stage_eta_poly_realgas = None
        if stage_exit is not None:
            stage_eta_poly = ((gamma - 1.0) / gamma) * math.log(stage_PR) / math.log(T02 / T01)
            # REAL-GAS exponent, with cp mass-averaged over the ACTUAL compression.
            # eta_poly is DEFINED through (gamma-1)/gamma, so an error in gamma lands in eta
            # ONE-FOR-ONE, with no physics in between. Air's cp rises 1005 -> 1026 J/kgK from
            # 288 to 484 K; mass-averaged over HECC's compression it is 1012.9, giving
            # gamma = 1.3954, not 1.3997. The constant-cp value above is therefore INFLATED
            # by ~0.8 eta points -- roughly HALF the efficiency deficit that was being hunted
            # for in the loss models. THIS is the number comparable to NASA's measurement.
            g_bar = fluid.gamma_mean(T01, T02)
            stage_eta_poly_realgas = (
                ((g_bar - 1.0) / g_bar) * math.log(stage_PR) / math.log(T02 / T01)
            )

        return OperatingPoint(
            PR=PR,
            psi=psi,
            stations=StationSeries(le_state, te_state),
            losses=losses_result,
            work_euler=work_euler,
            work_actual=work_actual,
            eta_poly=eta_poly,
            impeller_eta_is=impeller_eta_is,
            stage_states=tuple(stage_states),
            stage_PR=stage_PR,
            stage_eta_poly=stage_eta_poly,
            stage_eta_poly_realgas=stage_eta_poly_realgas,
        )
