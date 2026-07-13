"""Impeller loss models -- internal vs parasitic (docs/PHYSICS-RULES.md rule 6, defect D).

Upstream ships 16 centrifugal loss correlations in ``turbodesign/loss/compressor/otac.py``
that all declare ``LossType.Enthalpy`` -- a currency ``compressor_math.py`` silently drops
to ``Yp = 0``. The equations in those classes are frequently good physics; the wiring is
what's broken, and upstream's own translations carry bugs on top of that (its
``ImpellerDiscFrictionDaily`` divides by an inlet relative-velocity magnitude instead of
``mdot``; its ``ImpellerSkinFrictionJansen`` drops the ``L/d_h`` term Jansen's own equation
requires). This module does **not** import from ``otac.py`` (it depends on upstream
``BladeRow``, which this package deliberately never imports -- ``docs/PHYSICS-RULES.md``) and does not
reuse its arithmetic -- it re-derives each correlation from the primary/verified secondary
source in ``research/notes/02-loss-models-catalog.md``, cited inline below.

THE CENTRAL PHYSICS (docs/PHYSICS-RULES.md rule 6, research notes S0.1 -- Oh, Yoon & Chung 1997,
Kovář et al. 2021 Eq. 6)
------------------------------------------------------------------------------------
    INTERNAL  (incidence, skin friction, clearance, mixing, shock, choke):
        entropy generated INSIDE the blade passage. Destroys P0. Does NOT change the
        work input:                                   T02 UNCHANGED, P02 DOWN.
    PARASITIC (disc friction, recirculation, leakage):
        EXTRA SHAFT WORK that never becomes useful pressure rise. ADDS to the actual
        stagnation enthalpy rise, with no static-pressure benefit:
                                                        T02 UP, P02 UNCHANGED.

    eta = (dh_Euler - dh_internal - dh_static) / (dh_Euler + dh_parasitic)
                                                   ^^^^^^^^^^^^^^^^^^^^^^^^
                                     parasitic in the DENOMINATOR (it adds work)

Lumping the two into one coefficient (upstream's single ``Yp``) makes the efficiency
rollover unreproducible at any coefficient value -- see the module docstrings of
``tests/unit/test_slice3_internal_loss.py`` / ``test_slice4_parasitic_loss.py``.

HOW AN ENTHALPY LOSS BECOMES A PRESSURE LOSS (internal only)
--------------------------------------------------------------
Research notes S0.2: "Convert [Δh] to a stagnation-pressure loss via the isentropic
relation at the local state." The standard derivation (Denton, "Loss Mechanisms in
Turbomachines", ASME 93-GT-435, 1993 S2): for a loss occurring at ~constant static
enthalpy, the entropy generated satisfies ``T*ds ~= dh_loss``, so ``ds = dh_loss/T``
(evaluated at the local, i.e. impeller-exit, static temperature -- these are all
impeller-passage correlations). Internal losses do not change T0 (this is their
defining property), so at fixed T0R the ideal-gas entropy relation
``s = cp*ln(T/Tref) - R*ln(P/Pref)`` gives ``ds = -R*ln(P0R_actual/P0R_ideal)``, i.e.

    P0R2 = P0R2_is * exp(-ds_internal / R),    ds_internal = dh_internal_total / T2

``Stage.solve`` applies this INSIDE the impeller-exit continuity solve's residual, so
every trial Cm2 the root-finder tries recomputes the internal loss and the resulting
rho2 -- the reduced density from the pressure loss is genuinely COUPLED to the solved
Cm2, not a two-pass (solve-then-correct) approximation. That coupling is real: lower
rho2 -> higher Cm2 (same mdot) -> (via ``Vt2 = sigma*U2 + Cm2*tan(beta2b)``,
beta2b < 0) lower Vt2 -> lower Euler work -> T02 genuinely falls, by ~0.7 K at the HECC
design point. This is exactly the mechanism ``test_entropy_rises_across_an_internal_loss``
and ``test_internal_loss_lowers_pressure_ratio_and_efficiency`` are checking. What
internal loss still never does is add an explicit work TERM -- ``work_actual ==
work_euler`` bit-for-bit whenever no parasitic loss is configured
(``test_internal_loss_adds_no_work_term``); the T02 shift is entirely mediated through
Vt2, never a direct enthalpy addition.

Parasitic losses do not enter this relation at all (research notes S0.1: they raise
h02 directly, "with no static-pressure benefit") -- ``Stage.solve`` adds
``parasitic_total/cp`` to the STATION's reported T0 only, after the aerodynamic
(pressure/density/velocity-triangle) solve is complete, so P02 is bit-for-bit
independent of which parasitic models are configured.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import ClassVar, Dict, Literal, Optional, Protocol, Sequence, runtime_checkable

from .coefficients import FROZEN
from .state import Air


@dataclass(frozen=True)
class ImpellerLossState:
    """Everything an impeller loss correlation needs, at one (trial or converged)
    impeller-exit operating state. SI throughout: J/kg, m, m/s, K, Pa, kg/s, rad/s.

    ``rho2`` has two different meanings depending on WHEN a loss is evaluated
    (``turbodesign.centrifugal.solver.Stage.solve``'s ``_internal_losses_at``):

    - For **internal** losses (evaluated on every trial Cm2 inside the coupled
      continuity solve) it is the trial-state ISENTROPIC (loss-free) exit density,
      i.e. the density implied by ``P0R2_is`` at that Cm2 -- NOT the fully-converged,
      post-internal-loss density. This exists so a loss that itself needs a rho2/rho1
      compressibility ratio (Jansen clearance, S5) does not create a circular
      dependency on its own output; the gap between this and the true converged rho2
      is the ~2% internal-loss pressure correction documented in ``solver.py``,
      immaterial to a ratio living inside a square root.
    - For **parasitic** losses (evaluated once, after the coupled solve has
      converged) it is the exact converged exit density.

    It is ``None`` only when neither applies (e.g. a bare :class:`ImpellerLossState`
    built outside ``Stage.solve``, as in a unit test) -- a loss model that needs it
    must say so by raising (``AssertionError``/``TypeError``) rather than silently
    computing on ``None``.
    """

    fluid: Air
    mdot: float
    omega: float
    Z: float  # effective blade count, Impeller.Z_eff (splitter-aware). AN EXIT/PASSAGE
    # QUANTITY. Do NOT use it for anything evaluated at the INLET -- see Z_le.
    backsweep_deg: float  # beta2b, exit BLADE angle from meridional, negative = backswept

    # impeller inlet (station 1 / inducer eye): RMS streamline plus hub/shroud radii
    r1h: float
    r1s: float
    U1: float  # blade speed at the RMS streamline
    Vm1: float  # meridional velocity (assumed uniform across the span -- see delta_h docs)
    T1: float  # static temperature
    rho1: float  # static density

    # impeller exit (station 2)
    r2: float
    b2: float  # passage width along the exit cut, m
    tip_clearance: float  # Impeller.tip_clearance, eps, m -- shroud blade-tip gap
    U2: float
    Cm2: float
    Vt2: float
    W2: float
    T2: float  # static temperature
    rho2: Optional[float] = None  # None until the internal-loss pressure solve resolves it

    # Inducer blade metal angle at the RMS streamline, deg from meridional, and the LE
    # tangential blade thickness (Impeller.inducer_blade_angle_deg / .le_blade_thickness).
    # ``beta1b_deg`` is None when the caller's Impeller does not carry one --
    # ImpellerIncidenceConrad then RAISES rather than silently returning zero, which is
    # exactly the failure mode it used to have (see its docstring).
    beta1b_deg: Optional[float] = None
    t1: float = 0.0
    # BLADE COUNT AT THE LEADING EDGE -- the MAIN blades only (Impeller.n_blades), NOT
    # Z_eff.
    #
    # A splitter does not reach the inducer eye: HECC's splitter LE sits at r = 0.0675 m,
    # while the main-blade LE cut is at r_hub = 0.0380 m. At the station where inlet
    # blockage is evaluated, THE SPLITTERS ARE NOT THERE. Blocking the inlet circumference
    # with Z_eff = 25.4 instead of Z = 15 overstates the blockage by 69% (21.1% vs 12.3%),
    # which inflates Conrad's beta_opt from 49.19 to 52.17 deg and drives the design
    # incidence to +0.64 deg -- a spuriously "shockless" impeller, and a 30x error in the
    # incidence loss (2 J/kg where the geometry gives ~65).
    #
    # That is not a cosmetic difference: it is the difference between the incidence
    # minimum landing ON the design point (which would look like a triumphant validation)
    # and landing ABOVE it (which is the open anomaly this model actually has). A bug that
    # manufactures agreement is the most dangerous kind.
    #
    # Defaults to Z (the old behaviour) ONLY so that a caller who never set it keeps
    # working; every Impeller-driven path sets it explicitly in Stage.solve.
    Z_le: Optional[float] = None
    # Inducer throat area, m^2 (Impeller.throat_area). None -> ImpellerChokeAungier falls
    # back to the A_th := A1 stand-in, which is inert (debt D5).
    A_th: Optional[float] = None

    # ---- slice S1 (docs/centrifugal/17-tdd-plan.md) -- the camberline-bug fix ----
    #
    # TWO DIFFERENT length quantities, deliberately kept as two fields so nothing can
    # conflate them again (that conflation -- one closure feeding three different
    # consumers that each wanted something else -- is exactly what slice S1 fixed):
    #
    # L_blade (Impeller.l_blade_m): the through-blade / camberline length. Consumed by
    # ImpellerSkinFrictionJansen (the wetted-passage length Jansen's Cf acts over) and
    # ImpellerMixingAungier (the blade-to-blade Delta-W integrates over the same
    # length). None -> _camberline_length() falls back to the CORRECTED closure
    # (Gambini & Vellini's axial-length estimator, global flow coefficient -- see that
    # function's docstring), marked # UNVERIFIED (equation number unconfirmed, primary
    # paywalled). The SHIPPED closure's flow coefficient (phi = Vm1/U1, 12.5x too large)
    # is UNREACHABLE from this field: it was corrected in place, not routed around.
    #
    # L_meridional (Impeller.l_main_m): the MERIDIONAL length, LE->TE in the x-r plane
    # -- A DIFFERENT QUANTITY, always <= L_blade (a twisted blade is longer than its
    # meridional projection). Consumed ONLY by ImpellerLeakageAungier's Lm (research
    # notes 02 S10a): the leakage jet's driving pressure difference integrates along
    # the flow path, not along the twisted blade. Before this slice, leakage used
    # _camberline_length() (i.e. L_blade) as a PROXY for this -- "Lm := L_tilde", an
    # unsourced closure doing double duty (docs/centrifugal/12-model-as-implemented.md
    # S7.2). None -> ImpellerLeakageAungier falls back to _camberline_length(), the
    # SAME proxy as before, for a caller that has not supplied either measured length
    # (matching the Z_le/A_th "old behaviour by default" pattern above).
    L_blade: Optional[float] = None
    L_meridional: Optional[float] = None

    # ---- slice S2 (docs/centrifugal/17-tdd-plan.md) -- the Galvas K_BL/Z fix ----
    #
    # Galvas, NASA TN D-7487 (1973), Eq. (B59) + p. 5-6: D_f's own K_BL is 0.75 for a
    # "conventional" impeller and 0.6 "for impellers with splitters" (FORTRAN listing:
    # ``CONST1=0.75; IF(SPLT.EQ.1) CONST1=0.6``), and his Z is an INTEGER blade count
    # ("number of impeller blades at exit, Z3") -- never a splitter-weighted fraction.
    # The shipped code used the unconditional 0.75 together with ``Z = Z_eff``
    # (Aungier's FRACTIONAL, wetted-area-sense count) -- a hybrid that appears in no
    # source (data/coefficients.md S3.1, "D_f prefactor" row).
    #
    # Both fields default to the OLD (pre-slice-S2) behaviour when not supplied, exactly
    # like Z_le/A_th/L_blade above: a caller that never sets them keeps getting
    # K_BL = 0.75 and Z = state.Z (fractional Z_eff) -- see :func:`_diffusion_factor`.
    # ``Stage.solve`` sets both explicitly at BOTH ``ImpellerLossState`` construction
    # sites (the internal-loss trial state and the parasitic-block state) -- missing
    # one is the exact bug class this project keeps finding (Z_le, slice before this
    # one).
    #
    # Z_exit: the INTEGER blade count actually present at the impeller EXIT
    # (Impeller.n_blades + Impeller.n_splitters -- both blade rows reach the TE).
    # None -> _diffusion_factor falls back to state.Z (Z_eff, fractional).
    Z_exit: Optional[float] = None
    # has_splitters: selects Galvas's K_BL branch (0.6 if True, 0.75 if False/None).
    # For an UNSPLITTERED impeller (n_splitters == 0) this stays False and Z_exit
    # equals n_blades, which already equals Z_eff -- so D_f is BIT-IDENTICAL to the
    # pre-slice-S2 code on Eckardt rotors O/A (test_unsplittered_impeller_is_bit_identical,
    # tests/unit/test_diffusion_factor_provenance.py). That bit-identity matters: Kovář et al.
    # et al. 2021 validated this loss set against Eckardt with K_BL = 0.75, and this
    # slice must not silently change the physics being arbitrated against it.
    has_splitters: bool = False

    # ---- slice S4 (docs/centrifugal/17-tdd-plan.md) -- the exit metal blockage ----
    #
    # blockage2: Impeller.blockage, the SAME fraction the velocity-triangle solve
    # already applies to the exit area (``te_geom.with_blockage(...)`` in
    # ``Stage.solve``; ``docs/PHYSICS-RULES.md`` rule 4, ``area = 2*pi*r*b*(1-blockage)``).
    #
    # Before this slice, ``ImpellerLossState`` did not carry it at all:
    # :class:`ImpellerMixingAungier` computed its own ``A2 = 2*pi*r2*b2`` -- the
    # UNBLOCKED geometric area -- hardcoded, independent of whatever blockage the
    # triangle itself used. The moment ``Impeller.blockage != 0``, that is a silent
    # decoupling: the velocity triangle sees a BLOCKED exit area, the mixing loss sees
    # an UNBLOCKED one, and the two disagree about the geometry of the SAME station.
    # This field closes that gap by giving the loss model the one number it was
    # missing; :func:`ImpellerMixingAungier.delta_h` uses
    # ``A2 = 2*pi*r2*b2*(1-blockage2)`` instead.
    #
    # THIS IS GEOMETRIC (metal) blockage, not the aerodynamic (boundary-layer/wake)
    # blockage this project already killed once for being a circular free knob
    # (``Impeller.blockage``'s own docstring; ``data/coefficients.md``) -- see that
    # field's docstring and ``data/coefficients.md`` for the citations (Oh, Yoon &
    # Chung 1997 carries no blockage anywhere; Eckardt TM-75232's wake is
    # momentum-RICH, so a one-zone aerodynamic blockage over-corrects work by
    # construction). That door stays CLOSED; this field only ever carries the
    # GEOMETRIC (blade-metal) fraction ``Impeller.blockage`` was already given.
    #
    # Defaults to 0.0 -- the pre-slice-S4 behaviour -- so a caller that does not set
    # it (or an Impeller with ``blockage=0.0``, the field default) is BIT-IDENTICAL:
    # ``A2 = 2*pi*r2*b2*(1-0.0) == 2*pi*r2*b2`` exactly
    # (test_zero_blockage_is_bit_identical, tests/unit/test_exit_blockage.py).
    # ``Stage.solve`` sets it explicitly at BOTH ``ImpellerLossState`` construction
    # sites (the internal-loss trial state and the parasitic-block state) -- missing
    # one is the exact bug class this project keeps finding (Z_le, then Z_exit/
    # has_splitters, now this) --
    # test_both_loss_state_construction_sites_carry_blockage guards it directly.
    blockage2: float = 0.0

    @property
    def Z_inlet(self) -> float:
        """Blade count present AT THE LE. Main blades only -- splitters start downstream."""
        return self.Z if self.Z_le is None else self.Z_le

    @property
    def U1h(self) -> float:
        """Hub blade speed at the LE, omega*r1h -- exact, from geometry+omega."""
        return self.omega * self.r1h

    @property
    def U1s(self) -> float:
        """Shroud blade speed at the LE, omega*r1s -- exact, from geometry+omega."""
        return self.omega * self.r1s

    @property
    def W1h(self) -> float:
        """Hub relative velocity at the LE, hypot(Vm1, U1h).

        Assumes Vm1 (meridional velocity) is uniform across the inlet span -- the
        standard single-streamline meanline simplification (this module carries no
        hub/shroud meridional-velocity distribution, unlike Aungier's curvature-driven
        split in research notes S1c). Also assumes no inlet prewhirl (Vt1=0), matching
        every station this package's ``Stage.solve`` computes (slices 1-4 carry no IGV).
        """
        return math.hypot(self.Vm1, self.U1h)

    @property
    def W1s(self) -> float:
        """Shroud relative velocity at the LE -- see :attr:`W1h`."""
        return math.hypot(self.Vm1, self.U1s)


@runtime_checkable
class LossModel(Protocol):
    """The loss-model plugin seam. ``kind`` has NO default anywhere in this module --
    every concrete class below states it explicitly, because guessing it wrong
    inverts the physics (docs/PHYSICS-RULES.md rule 6).
    """

    kind: Literal["internal", "parasitic"]

    def delta_h(self, state: ImpellerLossState) -> float:
        """Specific enthalpy loss, J/kg (research notes S0.2's "enthalpy-loss form"),
        always >= 0 -- the sign of its effect on P0/T0 is determined by ``kind``, not
        by the sign of this return value.
        """
        ...


@dataclass(frozen=True)
class EvaluatedLosses:
    """The result of running a :class:`LossSet` against one operating state --
    resolved per-model contributions, keyed by class name, split by ``kind``.
    """

    internal: Dict[str, float]
    parasitic: Dict[str, float]

    @property
    def internal_total(self) -> float:
        return sum(self.internal.values())

    @property
    def parasitic_total(self) -> float:
        return sum(self.parasitic.values())


NO_LOSSES = EvaluatedLosses(internal={}, parasitic={})


def head_loss_correction(state: ImpellerLossState) -> float:
    """Aungier (1995) Eq. (32) -- ``f_c``, the head-loss (compressibility) correction.

        f_c = [2*T'_T2 / (T'_T1 + T'_T2)] * [2*(P'_T1 - P1) / (rho1*W1**2)]

    ADOPTED E-R11 (C3), PI sign-off (docs/centrifugal/41-prereg-adoption.md S5).
    Result: docs/centrifugal/42-r11-result.md.

    WHY IT EXISTS -- Aungier, p.363, verbatim:

        "the impeller loss models have been developed as adiabatic head (enthalpy)
         losses. This is strictly valid for incompressible flow, only. ... the losses
         occur in the blade passage, but are imposed at the impeller tip. ... The
         consequence of these two factors is a DETERIORATION IN PREDICTION ACCURACY AS
         ROTATIONAL MACH NUMBER (OR PRESSURE RATIO) INCREASES."

    and:

        "With this head loss correction procedure, the analysis has been used for stage
         pressure ratios up to 3.5 with no observable difference in prediction accuracy
         for different pressure ratios. WITHOUT this procedure, the deterioration of
         prediction accuracy with pressure ratio was quite noticeable."

    Both factors exceed 1 in a compressor and BOTH GROW WITH MACH/PR, so f_c raises the
    internal loss and LOWERS predicted PR, by more at higher PR.

    *** IT ADDS NO FREE PARAMETER. *** f_c is computed from thermodynamic state alone --
    there is nothing here to tune, so adopting it cannot be a fit. Measured AT THE DESIGN
    POINT: f_c = 1.1198 (Eckardt) / 1.2883 (HECC) / 1.2859 (CC3).

    ⚠️ THOSE ARE DESIGN-POINT VALUES ONLY. f_c VARIES WITH OPERATING POINT -- on HECC it
    falls monotonically with flow: 1.2883 (100 %) -> 1.2753 (90 %) -> 1.2648 (80 %) ->
    1.2562 (70 %). Never quote one number per machine as if it were a constant.

    ⚠️ AND IT IS NOT A POST-HOC MULTIPLIER. That is how it is APPLIED, not how it BEHAVES:
    more internal loss -> lower rho2 -> continuity must find a HIGHER Cm2 -> W_bar changes
    -> the correlations return DIFFERENT RAW NUMBERS before f_c ever multiplies them. HECC
    raw skin friction: 913.6 (control) -> 897.9 (f_c off) -> 912.1 (f_c on) J/kg -- it moved
    TWICE, IN OPPOSITE DIRECTIONS. Anyone reasoning about f_c as "just a x1.29 on the loss
    sum" is WRONG for every quantity downstream of rho2. See 42-r11-result.md S8.

    >>> THE ASSUMPTION, STATED WHERE IT IS MADE. <<<
    Aungier's Eq. (33) applies f_c to "the various loss coefficients PRESENTED IN THIS
    PAPER". OUR SET IS A HYBRID (Oh + Aungier + Jansen + Coppage + Conrad) that NO
    PUBLISHED PAPER EVER VALIDATED AS A SET. Applying f_c to it rests on the argument
    that Aungier's critique is of the FORMULATION -- adiabatic head losses imposed at the
    tip -- and not of his particular correlations; Oh's terms share that formulation
    exactly. THAT ARGUMENT IS OURS, NOT AUNGIER'S. It is an assumption, and it is
    reported as one (data/coefficients.md; docs/ASSUMPTIONS-AND-CORRELATIONS.md).

    AND IT DOES NOT RESCUE THE MODEL: it recovers 2.13 pp of HECC's 5.30 pp error and
    makes ECKARDT WORSE (-4.51 -> -5.14 %). Pre-registered as C3c and CONFIRMED. It is
    a MISSING TERM, not a unifying mechanism -- do not present it as one.
    """
    W1 = math.hypot(state.Vm1, state.U1)
    cp = state.fluid.cp_at(state.T1)

    TT1r = state.T1 + W1**2 / (2.0 * cp)  # relative total T, station 1
    TT2r = state.T2 + state.W2**2 / (2.0 * cp)  # relative total T, station 2

    P1 = state.rho1 * state.fluid.R * state.T1  # static P, station 1
    g = state.fluid.gamma_mean(state.T1, TT1r)
    PT1r = P1 * (TT1r / state.T1) ** (g / (g - 1.0))  # relative total P, station 1

    return (2.0 * TT2r / (TT1r + TT2r)) * (2.0 * (PT1r - P1) / (state.rho1 * W1**2))


@dataclass(frozen=True)
class LossSet:
    """A configured collection of loss models, e.g. ``LossSet([ImpellerSkinFrictionJansen()])``.

    ``Stage.solve`` (``turbodesign.centrifugal.solver``) evaluates this against the
    impeller's converged operating state and reports the result as ``OperatingPoint.losses``.

    ``apply_head_loss_correction`` switches Aungier's ``f_c`` (Eq. 32) on. It scales the
    INTERNAL sum ONLY -- Eq. (33) subtracts ``f_c*SUM(dq)`` from the BLADE work ``I_B``,
    while the parasitic terms are SEPARATE work inputs in his Eq. (1). Scaling parasitic
    work by f_c would ADD WORK that no source puts there, violating rule 6.

    It defaults to ``False`` here, and ``OhLossSet`` overrides it to ``True`` (E-R11).

    🛑 **THIS DOCSTRING USED TO CLAIM THAT** ``OhLossSet(apply_head_loss_correction=False)``
    **"is bit-identical to the set every pre-E-R11 number was produced with". THAT WAS
    FALSE**, and it was false in the most dangerous way available: it was TRUE WHEN WRITTEN.

    E-R11 adopted three corrections. C3 (``f_c``) is this flag -- so turning it off does undo
    C3. But **C1** (mixing: the EXIT relative tangential, Aungier Eq. 26) and **C2** (leakage:
    the missing ``* U2``) were merged **INTO** :class:`ImpellerMixingAungier` and
    :class:`ImpellerLeakageAungier` themselves. They are not behind a flag. They cannot be
    turned off. So ``OhLossSet(apply_head_loss_correction=False)`` reproduces **C1+C2**, not
    the control -- and anyone who trusted this docstring to regenerate the published control
    would have got the wrong model, with no error and nothing red.

    **The pre-adoption control is reconstructed in** ``my_scripts/preadoption_control.py``
    (the two closures recovered from ``af7088c``), and
    ``tests/validation/test_preadoption_control.py`` proves that reconstruction reproduces
    the published table cell for cell. **Use that. Do not use this flag alone.**
    """

    models: Sequence[LossModel]
    apply_head_loss_correction: bool = False

    def evaluate(self, state: ImpellerLossState) -> EvaluatedLosses:
        """Run every configured model and split the result by ``kind``.

        Raises ``ValueError`` on any ``kind`` other than "internal"/"parasitic" --
        a configured loss model that contributes to NEITHER bucket is exactly defect D
        (``test_a_configured_loss_model_cannot_silently_contribute_nothing``): silence
        is not an option here either.
        """
        internal: Dict[str, float] = {}
        parasitic: Dict[str, float] = {}
        for model in self.models:
            name = type(model).__name__
            dh = model.delta_h(state)
            if model.kind == "internal":
                internal[name] = dh
            elif model.kind == "parasitic":
                parasitic[name] = dh
            else:
                raise ValueError(
                    f"{name}.kind={model.kind!r} is neither 'internal' nor 'parasitic' -- "
                    "LossModel.kind has no default because guessing wrong inverts the "
                    "physics (docs/PHYSICS-RULES.md rule 6)."
                )

        if self.apply_head_loss_correction and internal:
            # INTERNAL ONLY -- see the class docstring and rule 6.
            fc = head_loss_correction(state)
            internal = {k: fc * v for k, v in internal.items()}

        return EvaluatedLosses(internal=internal, parasitic=parasitic)


# --------------------------------------------------------------------------- internal


def _global_flow_coefficient(state: ImpellerLossState) -> float:
    """The GLOBAL flow coefficient phi_t used by Gambini & Vellini's axial-length
    closure (see :func:`_axial_length_Lz`) -- docs/centrifugal/17-tdd-plan.md slice S1.

    The literature-review identification of the source's third term is
    ``1.58*(dt**2 - dh**2)*(C1m/U2)``. Read literally with ``dt``/``dh`` as DIAMETERS
    this is dimensionally wrong (m^2, not dimensionless) to sit inside a sum with
    ``0.014`` and ``0.023*r2/r1h``. It is dimensionally consistent, and numerically
    matches the reviewed value (0.0608 for HECC) to 3 s.f., if ``dt``/``dh`` are the
    NORMALIZED diameter ratios ``D1s/D2``, ``D1h/D2`` common in this literature
    (Whitfield & Baines-style notation) -- i.e.

        phi_t = ((D1s/D2)^2 - (D1h/D2)^2) * (C1m/U2)
              = ((r1s/r2)^2 - (r1h/r2)^2) * (Vm1/U2)

    verified numerically against the HECC design point: 0.06076 vs the reviewed 0.0608.
    Equation number **UNVERIFIED** (Gambini & Vellini, Springer 2021 Ch. 6, paywalled;
    see :func:`_axial_length_Lz`'s docstring) -- this function fixes the ARITHMETIC
    (which flow coefficient), not the citation.
    """
    r1h, r1s, r2 = state.r1h, state.r1s, state.r2
    return ((r1s / r2) ** 2 - (r1h / r2) ** 2) * (state.Vm1 / state.U2)


def _axial_length_Lz(state: ImpellerLossState) -> float:
    """L_z, the impeller axial-length estimate inside Gambini & Vellini's camberline
    closure (docs/centrifugal/17-tdd-plan.md slice S1; lineage debt D4, slice S0).

    ⚠️ **THE BUG, NOW FIXED.** This closure's third term is the GLOBAL flow coefficient
    (:func:`_global_flow_coefficient`, HECC: 0.0608), not a LOCAL one. The code used to
    compute ``phi = Vm1/U1`` (the LE flow coefficient, HECC: 0.759) instead --
    **12.5x too large**, in the dominant term. That made ``L_z = 0.580 m``, where
    NASA's own geometry gives a true hub LE->TE axial extent of **0.162 m**: the code
    believed the impeller was LONGER than it is WIDE (``D2 = 0.432 m``). The
    ``phi = Vm1/U1`` form is UNREACHABLE from this function -- there is no code path
    left that returns it.

    Only reached as a FALLBACK, when the caller has not supplied a measured
    :attr:`ImpellerLossState.L_blade` (see :func:`_camberline_length`) -- this closure
    is an *estimator for when geometry is unknown*; NASA Appendix C gives HECC's
    geometry directly, and a proxy for a derivable quantity is a defect (the same
    reasoning that retired ``A_th := A1`` and the shaft-centerline ``Z_eff`` proxy).
    Equation number **UNVERIFIED** (Gambini & Vellini, *Turbomachinery: Fundamentals,
    Selection and Preliminary Design*, Springer (2021), Ch. 6 -- book PAYWALLED).
    """
    r1h, r2 = state.r1h, state.r2
    phi_t = _global_flow_coefficient(state)
    c0 = FROZEN["impeller.axial_length.gambini_c0"]
    c1 = FROZEN["impeller.axial_length.gambini_c1"]
    c2 = FROZEN["impeller.axial_length.gambini_c2"]
    return 2.0 * r2 * (c0 + c1 * r2 / r1h + c2 * phi_t)


def _camberline_length(state: ImpellerLossState) -> float:
    """Blade through-blade (camberline) length, ``L_blade``.

    docs/centrifugal/17-tdd-plan.md slice S1: **prefer measured geometry, fall back to
    the corrected closure.** ``state.L_blade`` (``Impeller.l_blade_m``) is NASA
    Appendix C's own along-blade length for HECC -- the MEASURED 3-D camberline arc,
    0.237875 m (``my_scripts/extract_hecc_blade_angles.py``, ``data/coefficients.md``;
    NOT the REJECTED LE-angle projection, 0.292394 m -- see
    :class:`~turbodesign.centrifugal.components.Impeller`'s ``l_blade_m`` docstring
    for why) -- used directly, with no closure involved, whenever the caller supplies
    it.

    Absent that, falls back to the closure (:func:`_axial_length_Lz`,
    :func:`_global_flow_coefficient`), marked **# UNVERIFIED** (equation number
    unconfirmed, Gambini & Vellini's book paywalled) -- corrected in place: the old
    ``phi = Vm1/U1`` form is gone, not routed around.

    Shared by every loss model that needs a blade-passage length scale
    (:class:`ImpellerSkinFrictionJansen`'s skin friction, :class:`ImpellerMixingAungier`'s
    blade-to-blade Delta-W). ``ImpellerSkinFrictionJansen`` used to carry its OWN
    inline copy of this same closure, "left unfactored so its slice-3/4 test coverage
    never has to re-prove this refactor" -- that rationale EXPIRES the moment BOTH
    sites need to change (slice S1 fixes the formula), so slice S1 unifies them: both
    now call this one function, and there is exactly one place left to get it wrong.

    ``ImpellerLeakageAungier`` does **NOT** call this function for its own length scale
    any more: it wants the MERIDIONAL length (``state.L_meridional`` /
    ``Impeller.l_main_m``), a DIFFERENT quantity (docs/centrifugal/
    12-model-as-implemented.md S7.2) -- see that class's docstring.
    """
    if state.L_blade is not None:
        return state.L_blade

    r1h, r1s, r2, b2 = state.r1h, state.r1s, state.r2, state.b2
    beta2b = math.radians(abs(state.backsweep_deg))
    Vm1 = state.Vm1

    beta1h = math.atan2(state.U1h, Vm1)  # shockless-entry stand-in, as elsewhere
    beta1s = math.atan2(state.U1s, Vm1)

    L_z = _axial_length_Lz(state)
    return (
        (math.pi / 4.0)
        * (2.0 * r2 - (r1h + r1s) - b2 + 2.0 * L_z)
        / (math.cos(beta1h) + math.cos(beta1s) + math.cos(beta2b))
    )


def _diffusion_factor(state: ImpellerLossState) -> float:
    """Coppage/Galvas impeller diffusion factor D_f (research notes 02 S4a; Galvas
    NASA TN D-7487 (1973) "Blade loading loss", Eq. (B59); Kovář et al. 2021 Eq. (26),
    verified). Drives BOTH the Coppage blade-loading loss (S4a) and Oh's own
    recirculation loss (S9b) -- research notes 01 S4.2's table caption: "Diffusion
    factor D_f (drives #4 and #9)". Computed once here, shared, rather than re-derived
    per class.

        D_f = 1 - W2/W1s + (K_BL*dh_Euler/U2^2) /
              [ (W1s/W2) * ( (Z/pi)*(1-D1s/D2) + 2*D1s/D2 ) ]

    ``dh_Euler = U2*Vt2`` (no inlet swirl anywhere in this module's velocity
    triangles, matching ``turbodesign.centrifugal.solver.Stage.solve``).

    ``K_BL`` and ``Z`` -- SOURCED, slice S2 (docs/centrifugal/17-tdd-plan.md;
    data/coefficients.md S3.1). Galvas, Eq. (B59) + p. 5-6, verbatim:

        "a value of 0.75 is used for K_BL for conventional impellers and a value of
         0.6 is used for impellers with splitters ... A parametric study of
         calculated diffusion factors with a variation in the number of blades
         indicated that changing the constant to 0.6 would compensate for the
         changing solidity near the exit."

    and his FORTRAN listing: ``CONST1=0.75`` ... ``IF(SPLT.EQ.1) CONST1=0.6``, read
    directly into ``DF=``. His own input list defines Z as "number of impeller blades
    at exit, Z3" -- an INTEGER count. He *considered* varying Z to represent
    splitters and REJECTED it in favour of the K_BL constant; nowhere in Galvas does
    a fractional blade count appear inside D_f.

    The code this replaces was a hybrid that appears in NO source: Coppage/Galvas's
    D_f + the NO-SPLITTER constant 0.75, unconditionally, + Aungier's fractional
    ``Z_eff`` (``state.Z``) from a DIFFERENT loss framework (wetted area / hydraulic
    diameter -- ``Impeller.Z_eff``'s own docstring). ``state.Z_exit``/
    ``state.has_splitters`` (``ImpellerLossState``, slice S2) select the sourced
    (Z, K_BL) pair; both default to the OLD values (``state.Z``, 0.75) so a caller
    that does not populate them is UNCHANGED -- see ``ImpellerLossState.Z_exit``'s
    docstring for why that makes an unsplittered impeller (Eckardt O/A) bit-identical.
    """
    W2, W1s = state.W2, state.W1s
    U2, Vt2 = state.U2, state.Vt2
    D2, D1s = 2.0 * state.r2, 2.0 * state.r1s
    Z = state.Z if state.Z_exit is None else state.Z_exit
    K_BL = (
        FROZEN["impeller.diffusion_factor.k_bl_splittered"]
        if state.has_splitters
        else FROZEN["impeller.diffusion_factor.k_bl_unsplittered"]
    )

    dh_euler = U2 * Vt2
    denom = (W1s / W2) * ((Z / math.pi) * (1.0 - D1s / D2) + 2.0 * D1s / D2)
    return 1.0 - W2 / W1s + (K_BL * dh_euler / U2**2) / denom


@dataclass(frozen=True)
class ImpellerSkinFrictionJansen:
    """Impeller skin-friction loss -- Oh, Yoon & Chung (1997), Table 6, p. 336, row
    "Skin friction loss" (read from a rendered page image of the primary). Oh attributes
    the correlation to Jansen (1967); **JANSEN (1967) ITSELF IS NOT HELD** -- our authority
    for the printed form below is OH'S TABLE 6, and nothing else. Do not cite Jansen for it.
    (The class keeps its ``...Jansen`` name only because renaming it is a public-API change
    orthogonal to this fix; the SOURCE is Oh.) INTERNAL: turbulent wall friction over the
    wetted passage area destroys total pressure without adding work.

        Δh_sf = 2 * Cf * (L̃/d_h) * W̄²
        W̄ = (V_1t + V_2 + W_1t + 2*W_1h + 3*W_2) / 8       <- Oh 1997, Table 6

    **W̄ WAS WRONG, AND THE ERROR WAS IN THE SOURCE WE COPIED.** This class previously
    computed ``W̄ = (2*W2 + W1s - W1h) / 4`` -- a byte-faithful transcription of Kovář
    et al. 2021, *Energies* 14(24) 8545, Eq. (33). We did not mis-transcribe it; **Kovář's
    Eq. (33) is itself broken, and it fails on its own terms: the weights sum to 2 over a
    denominator of 4.** Put W1s = W1h = W2 = W into it and this "mean" returns W/2. *A mean
    that does not reproduce a constant is not a mean.* At the HECC design point it returned
    **101.95 m/s -- BELOW the minimum (145.3 m/s) of the three velocities it averages** --
    and skin friction came out low by ~3.7x.

    The repair is NOT a patch of Kovář (a "+W1h" 3-term variant normalises correctly but has
    no primary behind it). It is **Oh's printed five-term form**: weights 1+1+1+2+3 = 8 = the
    denominator, all terms POSITIVE, so it reproduces a constant exactly. Frames, as Oh
    prints them: ``V_1t`` and ``V_2`` are **ABSOLUTE** velocities (inducer tip / impeller
    exit), ``W_1t``/``W_1h`` are **RELATIVE** at the inducer tip/hub, ``W_2`` is **RELATIVE**
    at exit. With no inlet prewhirl anywhere in this module, the absolute inlet velocity is
    purely meridional, so ``V_1t = Vm1`` exactly (uniform-Vm1 across the span --
    :attr:`ImpellerLossState.W1h`); ``V_2 = hypot(Cm2, Vt2)``.

    This change was compelled by the normalisation failure, not chosen for agreement: it
    makes HECC better, Eckardt WORSE, and leaves CC3 falsified
    (docs/centrifugal/45-jansen-skinfriction-fix.md).

    ``d_h`` (Jansen's hydraulic diameter, notes S5.0 / Kovář et al. Eq. 29) is the average of the
    exit and inlet hydraulic diameters of one blade passage:

        d_h = D2 * [ cos(b2)/(Z/pi + 2*D2*cos(b2)/b2)
                    + avg(D1s/D2, D1h/D2)*avg(cos(b1s),cos(b1h))
                      / (Z/pi + (D1s+D1h)/(D1s-D1h)*avg(cos(b1s),cos(b1h))) ]

    Blade angles ``b1h``, ``b1s`` (inlet, hub/shroud) are not carried by ``Impeller`` --
    this module's ``Stage`` has no incidence model, i.e. it assumes shockless entry
    (flow angle == blade angle at the LE), so the LE relative-flow angles
    ``atan(U1h/Vm1)``/``atan(U1s/Vm1)`` (uniform-Vm1 approximation, see
    :attr:`ImpellerLossState.W1h`) stand in for them exactly -- not an extra
    approximation on top of the model's own assumptions, the SAME one.

    ``L̃`` (blade through-blade / camberline length) is :func:`_camberline_length`:
    NASA's own measured geometry when supplied (``Impeller.l_blade_m``, the MEASURED
    3-D camberline arc, 0.237875 m for HECC -- NOT the REJECTED LE-angle projection,
    0.292394 m), falling back to the corrected Gambini & Vellini closure
    (``# UNVERIFIED`` equation number) otherwise -- docs/centrifugal/17-tdd-plan.md
    slice S1.

    **Unified, slice S1.** This class used to carry its OWN inline copy of the
    camberline-length closure, "left unfactored so its slice-3/4 test coverage never
    has to re-prove this refactor" -- bit-for-bit identical to
    :func:`_camberline_length`'s computation, by construction. That rationale expired
    the moment the closure itself needed fixing (the shipped ``phi = Vm1/U1`` was
    12.5x too large, in the dominant term of ``L_z`` -- see
    :func:`_axial_length_Lz`'s docstring): a bug living in two copies must be fixed in
    two places, so slice S1 unifies them into the one call below.

    ``Cf = 0.004`` is Galvas's fixed input (NASA TN D-7487 (1973), "Skin friction loss",
    verified verbatim) -- notes S5.0: "a perfectly defensible default for a preliminary
    code", chosen over the UNVERIFIED Blasius-type fit CIMdes uses because it carries an
    actual primary citation. Overridable via the constructor for a case with a known
    Moody/Colebrook roughness.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    Cf: float = FROZEN[
        "impeller.skin_friction.cf"
    ]  # Galvas TN D-7487 (1973), verified verbatim; see class docstring.

    def delta_h(self, state: ImpellerLossState) -> float:
        r1h, r1s, r2, b2 = state.r1h, state.r1s, state.r2, state.b2
        Z = state.Z
        beta2b = math.radians(abs(state.backsweep_deg))
        Vm1 = state.Vm1
        W1h, W1s, W2 = state.W1h, state.W1s, state.W2

        beta1h = math.atan2(state.U1h, Vm1)  # shockless-entry stand-in; see class docstring
        beta1s = math.atan2(state.U1s, Vm1)

        D2, D1h, D1s = 2.0 * r2, 2.0 * r1h, 2.0 * r1s
        cos_b2 = math.cos(beta2b)
        cos_b1_mean = 0.5 * (math.cos(beta1s) + math.cos(beta1h))

        # E-R3/D1 (docs/centrifugal/26-prereg-r3.md): this denominator previously carried a
        # spurious `2.0 *` on the D2 term. Kovar et al. 2021 Eq. (29) does NOT have it, and
        # `term2` below -- the inlet counterpart, same construction -- does not have it
        # either: the code was inconsistent with itself as well as with its declared source,
        # while the ledger row called the transcription "VERBATIM". Removing it raises d_h,
        # lowers skin friction, and moves HECC's PR error +5.79% -> +6.06%: FURTHER FROM
        # NASA. Reported, not tuned.
        # CAVEAT, against our own interest: Jansen (1967) itself is NOT HELD. If it is ever
        # acquired and DOES carry the 2, this change is reversed, and the reversal is reported.
        term1 = cos_b2 / (Z / math.pi + D2 * cos_b2 / b2)
        term2 = (0.5 * (D1s / D2 + D1h / D2) * cos_b1_mean) / (
            Z / math.pi + (D1s + D1h) / (D1s - D1h) * cos_b1_mean
        )
        d_h = D2 * (term1 + term2)

        L_tilde = _camberline_length(state)  # slice S1: unified, no more inline copy

        # Oh, Yoon & Chung (1997) Table 6: W_bar = (V_1t + V_2 + W_1t + 2*W_1h + 3*W_2)/8.
        # Weights 1+1+1+2+3 = 8 = the denominator, so W_bar == W when all five are W. The
        # form this replaces -- Kovar et al. 2021 Eq. (33), (2*W2 + W1s - W1h)/4 -- summed
        # its weights to 2 over a denominator of 4 and returned W/2 for a constant field.
        # FRAMES (Oh's symbols): V_* ABSOLUTE, W_* RELATIVE. No inlet prewhirl in this
        # module => the absolute inlet velocity is purely meridional, so V_1t = Vm1 exactly.
        V1t = Vm1  # ABSOLUTE at the inducer tip (Vt1 = 0; uniform-Vm1, see W1h's docstring)
        V2 = math.hypot(state.Cm2, state.Vt2)  # ABSOLUTE at the impeller exit
        W_bar = (V1t + V2 + W1s + 2.0 * W1h + 3.0 * W2) / 8.0

        dh = 2.0 * self.Cf * (L_tilde / d_h) * W_bar**2
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerIncidenceConrad:
    """Conrad (1979) inducer incidence loss -- Oh, Yoon & Chung (1997)'s selection
    (research notes 02 S1b; Kovář et al. 2021 Eq. (13), verified). INTERNAL: the
    relative-velocity component normal to the optimum-incidence direction is assumed
    destroyed.

        Δh_inc = f_inc * W*^2 / 2,      W* = W_1xi * sin|beta_opt - beta_1xi|
        beta_opt = arctan[ (pi*D_1xi / (pi*D_1xi - Z*t_1xi)) * tan(beta_tilde_1xi) ]

    (``beta_opt`` is the Galvas/Stanitz blockage-corrected optimum angle, research
    notes 02 S1a.) ``f_inc = 0.6`` is Kovář et al. 2021 Eq. (13)'s "commonly 0.6" mid-point
    of the published 0.5-0.7 range -- selected by Oh 1997 and by Kovář et al.'s best sets
    1759/2042/2137.

    **RESURRECTED (item 5e).** This model used to be structurally dead. With no inducer
    blade-metal-angle field on ``Impeller``, it had to define ``beta_tilde_1xi :=
    beta_1xi`` -- the blade angle set equal to the FLOW angle -- which makes
    ``beta_opt == beta_1xi`` identically, ``W* == 0``, and the loss EXACTLY ZERO at
    every operating point, forever. Its ``f_inc`` was a phantom parameter. That is the
    same silent-zero defect (``docs/PHYSICS-RULES.md`` defect D; upstream's ``LossType.Enthalpy ->
    Yp = 0``) that this module was written to avoid, reproduced inside it.

    The defence at the time was that a real inducer is twisted so the blade matches the
    flow at design, so a near-zero incidence loss AT a design point is expected. True --
    but it made the loss zero at EVERY point, including off-design, and incidence is
    precisely what makes efficiency roll over on a speedline. A model that cannot be
    wrong cannot be evidence.

    ``beta_1b`` is now supplied (``Impeller.inducer_blade_angle_deg``). It is not
    tabulated by NASA -- Table 2 gives only LEAN angles -- but it is DERIVABLE from
    NASA's own Appendix C blade coordinates: see ``my_scripts/extract_hecc_blade_angles.py``
    and ``data/hecc/blade_angles.csv`` (HECC, RMS streamline: 45.46 deg).

    If no blade angle is supplied this now RAISES. It does not fall back to the
    flow-angle stand-in, because that fallback is indistinguishable from "no loss" and
    a loss model that cannot affect the answer is worse than one that errors.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    f_inc: float = FROZEN[
        "impeller.incidence.f_inc"
    ]  # Kovář et al. 2021 Eq. (13); Oh 1997's / Kovář et al.'s best-sets choice

    def delta_h(self, state: ImpellerLossState) -> float:
        if state.beta1b_deg is None:
            raise ValueError(
                "ImpellerIncidenceConrad needs an inducer blade metal angle "
                "(Impeller.inducer_blade_angle_deg). Without it, beta_opt collapses to "
                "the flow angle, W* == 0, and this loss is identically zero at every "
                "operating point -- a dead model that the tests would certify as present "
                "(docs/PHYSICS-RULES.md defect D). Derive it from the blade coordinates: see "
                "my_scripts/extract_hecc_blade_angles.py."
            )
        Vm1, U1 = state.Vm1, state.U1
        beta_1xi = math.atan2(U1, Vm1)  # LE relative FLOW angle, from meridional
        beta_1b = math.radians(state.beta1b_deg)  # LE relative BLADE angle -- independent

        # Galvas/Stanitz blockage-corrected optimum-incidence angle. The blades block part
        # of the inlet circumference, so the flow that enters cleanly is turned slightly
        # more than the metal angle. With t1 = 0 this collapses to beta_opt == beta_1b --
        # the correct boundary case (NASA's tabulated sections meet at a point at the LE),
        # not a fudge.
        D_1xi = 2.0 * math.hypot(state.r1h, state.r1s) / math.sqrt(2.0)  # RMS diameter
        # Z_inlet, NOT Z. The splitters do not reach the inducer eye, so they cannot block
        # it (ImpellerLossState.Z_le). Using the splitter-aware Z_eff here overstated the
        # blockage by 69% and made the impeller look shockless at design.
        Z_le = state.Z_inlet
        blocked = math.pi * D_1xi - Z_le * state.t1
        if blocked <= 0.0:
            raise ValueError(
                f"inlet blade blockage exceeds the full circumference: Z_le*t1 = "
                f"{Z_le * state.t1:.4f} m >= pi*D_1xi = {math.pi * D_1xi:.4f} m"
            )
        beta_opt = math.atan((math.pi * D_1xi / blocked) * math.tan(beta_1b))

        W1xi = math.hypot(Vm1, U1)
        W_star = W1xi * math.sin(abs(beta_opt - beta_1xi))
        dh = self.f_inc * 0.5 * W_star**2
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerBladeLoadingCoppage:
    """Coppage (1956) blade-loading loss, via the diffusion factor -- Oh, Yoon &
    Chung (1997)'s selection (research notes 02 S4a; Galvas NASA TN D-7487 (1973)
    "Blade loading loss", verified verbatim: ``Δh_BL = 0.05 D_f^2 u_3^2``; ``D_f``
    from Kovář et al. 2021 Eq. (26), verified). INTERNAL: blade-to-blade pressure-gradient
    driven boundary-layer growth / secondary flow / separation.

        Δh_bl = 0.05 * D_f^2 * U2^2

    Coefficient: 0.05 (Galvas, verified verbatim). ``D_f`` from :func:`_diffusion_factor`.

    ==================== ⭐ THE LARGEST INTERNAL LOSS IN THE MODEL ====================
    Post-E-R11 (the mixing loss's correct station annihilates it), THIS TERM IS THE
    LARGEST INTERNAL LOSS ON ALL THREE MACHINES:

        Eckardt 1788 J/kg (58.2 %)   HECC 4295 (61.4 %)   CC3 3950 (66.8 %)

    *** AND ITS AUTHOR DISOWNS THE EXPRESSION IT COMES FROM. ***

    PRIMARY, HELD and read from the page image -- Coppage et al., WADC TR 55-257 (1956),
    p.11, S2.4, VERBATIM:

        "In the cases of the DIFFUSION AND BLADE LOADING LOSS, and the recirculation
         loss, EXPERIMENTAL RESULTS WERE TOO MEAGER TO BE OF MUCH ASSISTANCE. The
         expressions for these losses CAN ONLY BE REGARDED AS HYPOTHESES."

    And the coefficient itself, Eq. (2.34), p.16 -- note the lead-in:

        "The expression for Dq_DBL which is proposed, BASED ON AiRESEARCH EXPERIENCE, is"
            Dq_DBL = 0.050 * Delta**2                                          (2.34)

    NO data. NO figure. NO validity range. NO uncertainty. Anywhere.

    ⚠️ SCOPE -- DO NOT OVERQUOTE HIM. Coppage disowns THE EXPRESSIONS FOR TWO LOSSES
    (diffusion-and-blade-loading, and recirculation) -- NOT "his coefficients" as a class.
    In the very next paragraph he DEFENDS disk friction and skin friction as "already
    well-established". BLADE LOADING IS ONE OF THE TWO HE NAMES, so the claim STANDS --
    but it must be quoted this way and no other.

    ⚠️ NAMING: Coppage calls this the "DIFFUSION AND blade loading" loss (Dq_DBL) -- it
    lumps Lieblein-style diffusion together with blade loading in ONE term. This class's
    name drops half of that. The number is unaffected.

    See docs/ASSUMPTIONS-AND-CORRELATIONS.md (layer B) and data/coefficients.md (E5).
    ==================================================================================
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    k: float = FROZEN["impeller.blade_loading.k"]  # Galvas TN D-7487 (1973), verified verbatim

    def delta_h(self, state: ImpellerLossState) -> float:
        Df = _diffusion_factor(state)
        dh = self.k * Df**2 * state.U2**2
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerClearanceJansen:
    """Jansen (1967) impeller tip-clearance loss -- Oh, Yoon & Chung (1997)'s
    selection, and the default recommendation independent of Oh (research notes 02
    S6a; Kovář et al. 2021 Eq. (37), verified; confirmed independently by CIMdes
    ``loss.py::ClearanceLoss``, ``radcomp/impeller.py:98`` (credited "Jansen, Brasz"),
    and Yang, Liu & Zhao 2023 *Machines* 11(1):118 Eq. (8) (citation corrected, slice S0
    -- previously misattributed as "Li 2023"; same paper as the Z_eff formula below) --
    four independent confirmations). INTERNAL: flow
    leaking over the unshrouded blade tip from pressure to suction side.

        Δh_cl = 0.6 * (eps/b2) * |Cu2| * sqrt[
            (4*pi / (b2*Z)) * (r1s^2 - r1h^2) / ((r2-r1s)*(1+rho2/rho1)) * |Cu2| * Cm1xi
        ]

    Coefficient: 0.6 (verified, four independent sources). ``eps`` is
    ``Impeller.tip_clearance`` (this slice's new field). Guarded on ``|Cu2|`` -- the
    raw form goes complex at negative swirl (research notes 02 S6a).

    ``rho2`` here is the TRIAL ISENTROPIC exit density (see :class:`ImpellerLossState`'s
    docstring for why: this loss's own rho2/rho1 ratio would otherwise be circular
    with the internal-loss total it is part of). The gap to the fully-converged
    density is the ~2% internal-loss pressure correction documented in
    ``turbodesign.centrifugal.solver`` -- immaterial to a ratio inside a square root.

    ``Z`` -- **R9 (docs/centrifugal/20-review-fixes.md), the SEVENTH cross-station
    instance.** Jansen's own clearance formula specifies a blade COUNT (his eq. counts
    leakage paths around the shroud gap), not Aungier's fractional, wetted-length
    ``Z_eff``. At the impeller EXIT both blade rows (main + splitter) physically exist
    -- HECC: 30, not ``Z_eff = 25.4048``. ``state.Z_exit`` (already populated by
    ``Stage.solve`` for :func:`_diffusion_factor`'s slice-S2 fix) is read here too;
    ``None`` -> falls back to ``state.Z`` (``Z_eff``), so an unsplittered impeller
    (``n_splitters = 0`` => ``Z_exit == n_blades == Z_eff``) is bit-identical.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    jansen_coeff: float = FROZEN[
        "impeller.clearance.jansen_coeff"
    ]  # Jansen (1967); four confirmations

    def delta_h(self, state: ImpellerLossState) -> float:
        assert state.rho2 is not None, (
            "ImpellerClearanceJansen needs an exit-density estimate -- Stage.solve "
            "must supply the trial isentropic rho2 to internal-loss evaluation (see "
            "this class's and ImpellerLossState's docstrings)."
        )
        r1h, r1s, r2, b2 = state.r1h, state.r1s, state.r2, state.b2
        Z = state.Z if state.Z_exit is None else state.Z_exit
        Cu2 = state.Vt2
        Cm1xi = state.Vm1
        rho_ratio = state.rho2 / state.rho1
        eps = state.tip_clearance

        inner = (
            (4.0 * math.pi / (b2 * Z))
            * (r1s**2 - r1h**2)
            / ((r2 - r1s) * (1.0 + rho_ratio))
            * abs(Cu2)
            * Cm1xi
        )
        dh = self.jansen_coeff * (eps / b2) * abs(Cu2) * math.sqrt(max(inner, 0.0))
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerMixingAungier:
    """Aungier jet/wake mixing loss -- the single-zone-friendly alternative to
    Johnston & Dean (1966), used HERE INSTEAD OF Oh 1997's own literal choice
    (research notes 02 S7b vs S7a; Kovář et al. 2021 Eqs. (43)-(45), verified).

    ============================================================================
    ✅  THE STATION IS CORRECTED, AND THE PRIMARY IS NOW HELD. ADOPTED E-R11 (C1).
    ============================================================================
    Aungier (1995), ASME J. Turbomach. 117(3):360-366 -- HELD -- Eq. (26) prints:

        W_OUT**2 = [C_m2*A2/(pi*d2*b2)]**2 + W_u**2

    with a BARE ``W_u``: *** NO STATION SUBSCRIPT AT ALL. *** Not 1, not 2, not
    "1xi". Every OTHER symbol in the equation is at station 2 (C_m2, d2, b2, and
    A2 = "the discharge area inside the blades"), and W_SEP = W2. A station-1
    tangential inside an all-station-2 mixed-out state is incoherent.

    ==> Kovář et al. (2021) Eq. (45)'s ``W_u1xi`` is AN ADDITION KOVÁŘ MADE, NOT A
        TRANSCRIPTION. This code was byte-faithful to it for months, and its
        docstring called that transcription "verified" -- IT WAS VERIFIED AGAINST
        SOMETHING THE SOURCE NEVER SAID.

        BEING BYTE-FAITHFUL TO A PUBLISHED SOURCE IS NOT THE SAME AS BEING RIGHT.
        ONLY A SECOND SOURCE CAN CATCH A WRONG ONE.

    Four independent secondaries (Güllü 2015 METU; Sanz Solaesa 2016 KTH; Gong &
    Chen 2014; Harrison 2020 Purdue Table 2.10) all place it at station 2 -- and
    were right about the STATION while three of them were WRONG about the ALGEBRA
    (they print a meridional difference; Aungier prints a relative one). A
    head-count of secondaries got one axis right and the other backwards.
    See docs/centrifugal/40-aungier-primary-result.md.

    *** THIS CORRECTION ANNIHILATES THE TERM -- AND THAT IS THE RESULT. ***
    MEASURED at each machine's REAL blockage (E-R11, docs/centrifugal/42-r11-result.md S7):

        dh_mix:  Eckardt 2769 -> 2.15 J/kg    HECC 1735 -> 0.54    CC3 1198 -> 0.00

    Exit blockage does NOT rescue it: W_sep = hypot(Cm2, W_u2) vs
    W_out = hypot(Cm2*(1-B2), W_u2), and in a backswept impeller W_u2 >> Cm2, so
    trimming the MERIDIONAL leg by a few percent moves the HYPOTENUSE by a
    SECOND-ORDER amount -- which is then SQUARED.

    ==> AUNGIER'S MIXING LOSS, AT THE STATION HIS OWN EQUATION IMPLIES, IS
        STRUCTURALLY NEAR-ZERO FOR ANY ATTACHED IMPELLER. It fires ONLY through the
        D_eq > 2 STALL branch. IT IS A STALL TERM WEARING A MIXING TERM'S NAME --
        the ELEVENTH "term that cannot fire", and this one is IN THE CANON.

    And it makes the model WORSE: HECC's stage PR error goes +5.30 % -> +7.21 %.
    IT IS ADOPTED ANYWAY, because it is what the primary says. Adopting a change
    that moves the model AWAY from the measurement is the whole point of a control.

    Pre-registered: docs/centrifugal/41-prereg-adoption.md (locked, committed
    BEFORE the counterfactual). PI sign-off recorded there.
    H_doctrine SURVIVES: signs unchanged (-1.68 / +7.23 / +20.86 %).
    ============================================================================

    **Deviation from the literal Oh 1997 set, documented.** Johnston & Dean's
    ``Δh_mix`` needs a wake-width zeta and wake-mass-fraction lambda from Oh's own
    two-zone (jet/wake) split procedure (research notes 02 S7a: "requires a two-zone
    model, which is exactly why this loss is awkward in a single-zone code"). This
    module carries one streamline per station -- there is no second zone to split
    into, and no primary-sourced way to manufacture zeta/lambda from a single-zone
    state. Inventing them from an unsourced correlation would be exactly the kind of
    free knob this project has already had to kill twice (module docstring of
    ``tests/unit/test_slice5_full_loss_set.py``). The catalog's OWN guidance for this
    situation is explicit: "If you are single-zone, prefer 7b" (research notes 02
    S7a), and separately "the practical single-zone choice; all three of Kovář et al.'s best
    sets use it" (S7b). This class follows that guidance. INTERNAL: mixing out the
    jet/wake non-uniformity at the impeller exit.

        Δh_mix = 0.5 * (W_sep - W_out)^2
        W_sep = W2                if D_eq <= 2
              = W2 * D_eq / 2     if D_eq > 2
        D_eq = W_max / W2,   W_max = (W1xi + W2 + dW) / 2
        dW = 2*pi*D2*Cu2 / (Z * L_tilde)                        (blade-to-blade Delta-W)
        W_out = sqrt( (Cm2*A2 / (pi*D2*b2))^2 + Wu1xi^2 ),  Wu1xi = U1 (no inlet swirl)

    ``L_tilde`` is the THROUGH-BLADE (camberline) length, :func:`_camberline_length` --
    the SAME length :class:`ImpellerSkinFrictionJansen` uses (NASA's measured
    ``Impeller.l_blade_m`` when supplied, else the corrected closure, ``# UNVERIFIED``;
    docs/centrifugal/17-tdd-plan.md slice S1). Before slice S1 this closure's dominant
    term used the wrong (LE, not global) flow coefficient, 12.5x too large -- see
    :func:`_axial_length_Lz`.

    Coefficients: 0.5; the ``D_eq = 2`` (Lieblein) stall threshold (verified).

    ``A2`` -- **BLOCKED, slice S4 (docs/centrifugal/17-tdd-plan.md).** Before this
    slice this was hardcoded as the ideal (blockage-free) geometric area
    ``2*pi*r2*b2``, independent of whatever exit blockage the velocity-triangle solve
    itself used (``Stage.solve``'s ``te_geom.with_blockage(self.impeller.blockage)``).
    That is a silent decoupling the moment ``Impeller.blockage != 0``: the triangle
    computes ``Cm2`` against a BLOCKED area, this loss re-inflates it back to the
    UNBLOCKED one for ``W_out``, and the two disagree about the area of the SAME
    station. ``A2 = 2*pi*r2*b2*(1 - blockage2)`` (``ImpellerLossState.blockage2``,
    populated from ``Impeller.blockage`` at BOTH ``Stage.solve`` construction sites)
    closes that gap. At ``blockage2 = 0.0`` (the field default, and the pre-slice-S4
    behaviour) this is bit-for-bit the old ``2*pi*r2*b2``
    (test_zero_blockage_is_bit_identical, tests/unit/test_exit_blockage.py). At
    ``A2 = pi*D2*b2*(1-blockage2)`` the ``Cm2*A2/(pi*D2*b2)`` term reduces to
    ``Cm2*(1-blockage2)`` exactly, left unsimplified below so the formula stays
    traceable to the catalog.

    ``Z`` in ``dW`` -- **R9 (docs/centrifugal/20-review-fixes.md), the SEVENTH
    cross-station instance.** ``dW = 2*pi*D2*Cu2/(Z*L_tilde)`` is Aungier's
    blade-to-blade circulation term at the EXIT, where both blade rows (main +
    splitter) physically exist -- HECC: 30, not the fractional, wetted-length
    ``Z_eff = 25.4048``. ``state.Z_exit`` (already populated by ``Stage.solve`` for
    :func:`_diffusion_factor`'s slice-S2 fix) is read here too; ``None`` -> falls
    back to ``state.Z`` (``Z_eff``), so an unsplittered impeller
    (``n_splitters = 0`` => ``Z_exit == n_blades == Z_eff``) is bit-identical.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    mixing_k: float = FROZEN[
        "impeller.mixing.k"
    ]  # Aungier (2000) via Kovář et al. 2021 Eqs. (43)-(45)
    deq_threshold: float = FROZEN[
        "impeller.mixing.lieblein_deq_threshold"
    ]  # Lieblein stall threshold

    def delta_h(self, state: ImpellerLossState) -> float:
        r2, b2 = state.r2, state.b2
        Z = state.Z if state.Z_exit is None else state.Z_exit
        D2 = 2.0 * r2
        Cu2 = state.Vt2
        W2 = state.W2
        W1xi = math.hypot(state.Vm1, state.U1)
        L_tilde = _camberline_length(state)

        dW = 2.0 * math.pi * D2 * Cu2 / (Z * L_tilde)
        W_max = (W1xi + W2 + dW) / 2.0
        Deq = W_max / W2
        W_sep = W2 if Deq <= self.deq_threshold else W2 * Deq / self.deq_threshold

        # BLOCKED area, slice S4 -- see class docstring. blockage2 defaults to 0.0, so
        # this is bit-identical to the old `2*pi*r2*b2` whenever a caller (or an
        # Impeller with blockage=0.0) does not populate it.
        A2 = 2.0 * math.pi * r2 * b2 * (1.0 - state.blockage2)
        # THE EXIT relative tangential -- Aungier (1995) Eq. (26). ADOPTED E-R11 (C1).
        Wu2 = state.U2 - Cu2
        W_out = math.hypot(state.Cm2 * A2 / (math.pi * D2 * b2), Wu2)

        dh = self.mixing_k * (W_sep - W_out) ** 2
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerChokeAungier:
    """Aungier inducer choke loss -- Kovář et al. 2021 Eqs. (19)-(21) (verified
    transcription of Aungier 2000). NOT in Oh 1997's canonical set; required here
    because Oh's omission "largely overestimat[es] efficiency near choke point"
    (Kovář et al. 2021 S7.1, verified) and this project sweeps to choke in slice 10.
    INTERNAL: onset loss as the inducer throat approaches sonic.

        Δh_ch = W1xi * (0.05*x + x^7) / 2     if x > 0,  else 0
        x = 10 * (1.1 - Cr*A_th / A*_th),     Cr = sqrt(A1*sin(beta_1xi) / A_th)

    Coefficients: 0.05; the 7th power; onset scale 10; onset threshold 1.1; the
    SQUARE ROOT in ``Cr`` (verified, R7 (``docs/centrifugal/20-review-fixes.md``),
    against the held PDF ``research/papers/kovar_2021_energies-14-8545_loss-model-set.pdf``,
    Eq. (21): ``Cr = sqrt(A1 sin(beta~1xi) / Ath)`` -- the square root is present in
    the source and had been DROPPED here. See the note below the reciprocal-ratio
    discussion.

    ⚠ **Ratio inverted relative to the catalog's literal transcription, flagged --
    not a coefficient change.** Research notes 02 S2a prints
    ``x = 10*(1.1 - A*_th/(Cr*A_th))``. Implemented literally, that ratio is <1 for
    EVERY subsonic operating point (it -> 0 as flow -> 0, -> 1 only exactly at
    choke), so ``x > 0`` -- and the onset polynomial fires -- almost everywhere,
    including deep subsonic conditions (checked numerically at the HECC design point,
    M1s,rel = 0.84: the literal form gives ``x = +6.5``, ``Delta h_ch`` ~5.5e7 J/kg,
    an obvious blow-up). That contradicts the plain-English sentence next to the
    equation in the SAME source ("Only active when the inducer throat is near or
    past sonic") and the catalog's own separate flag on this equation ("check ...
    before trusting the magnitude" -- notes 02 S2a; the source PDF's equation here is
    exactly the kind of vector-rendered object notes 02 S14 warns transcribes
    unreliably). The reciprocal ratio, ``Cr*A_th/A*_th`` -- the AVAILABLE throat
    area over the REQUIRED sonic area, diverging to infinity at zero flow and
    dropping to 1 (x -> 1 > 0, loss switches on) exactly as the flow approaches
    choke -- matches the prose exactly and gives ``x = -11.1`` (loss = 0) at the
    HECC design point, which is correctly subsonic. This class uses the reciprocal
    form. Every published NUMBER (0.05, 10, 1.1, the 7th power) is unchanged.

    ⚠ **R7, docs/centrifugal/20-review-fixes.md: the square root in ``Cr`` had been
    DROPPED, under the claim above that "every published NUMBER is unchanged."**
    Kovář et al. 2021 Eq. (21), verified against the held PDF
    (``research/papers/kovar_2021_energies-14-8545_loss-model-set.pdf``): ``Cr = sqrt(A1 *
    sin(beta_1xi) / A_th)``. The code had ``Cr = A1*sin(beta_1xi)/A_th`` -- no
    ``sqrt`` -- which is not a coefficient error but a dropped operator. Restored.
    **Verified inert either way** (see the module-level sweep note referenced by
    ``data/coefficients.md``): with or without the square root, ``x < 0`` (loss =
    0) at every mass flow this model can reach on HECC (3.45-5.82 kg/s swept,
    including the model's OWN choke flow, ``inducer_choke_mdot`` -- the
    **diffuser** throat chokes first at 5.86 kg/s, before the inducer criterion
    ever crosses zero). A term that never fires still cannot be caught being
    wrong by this project's own test suite -- restoring the square root is
    therefore, and will remain, output-invariant on every case this repository
    currently validates against, not evidence the fix was unnecessary.

    ``A*_th`` (sonic/critical throat area) is the standard 1-D choking relation,
    evaluated at the SHROUD (research notes 01 S1: the shroud is "what chokes"),
    using the relative-frame stagnation state there and the absolute stagnation
    speed of sound at the LE RMS station -- matching the catalog's own mixed-frame
    form (``rho_0,rel a_0,rel`` prefactor, ``a_0`` -- absolute -- inside the bracket):

        A*_th = mdot / (rho0R_1s * a0R_1s) *
                [ (2 + (gamma-1)*(U1s/a0)^2) / (gamma+1) ] ** (-(gamma+1)/(2*(gamma-1)))

    **Approximation, flagged.** This module has no independent inducer-throat area
    (research notes 01 S1 calls "station 1th" "a sub-station, not always exposed").
    ``A1``, the LE annulus geometric area, is recovered EXACTLY from LE continuity
    (``A1 = mdot/(rho1*Vm1)``, since that is how Vm1 was solved for) and stands in
    for ``A_th`` -- which makes ``Cr = sqrt(sin(beta_1xi))`` exactly. The SAME throat-less
    simplification is used by :class:`ImpellerEntranceDiffusionAungier`. Away from
    choke (every operating point this slice tests) ``x < 0`` and the loss is zero by
    construction -- the expected result AT a design point, not evidence of a broken
    model. It becomes load-bearing, and needs a real throat station, once slice 10
    sweeps to choke.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    k_low: float = FROZEN[
        "impeller.choke.k_low"
    ]  # Aungier (2000) via Kovář et al. 2021 Eqs. (19)-(21)
    high_power: float = FROZEN["impeller.choke.high_power"]
    onset_scale: float = FROZEN["impeller.choke.onset_scale"]
    onset_threshold: float = FROZEN["impeller.choke.onset_threshold"]

    def delta_h(self, state: ImpellerLossState) -> float:
        fluid = state.fluid
        gamma, R, cp = fluid.gamma, fluid.R, fluid.cp
        Vm1, U1, U1s = state.Vm1, state.U1, state.U1s
        T1 = state.T1
        P1 = state.rho1 * R * T1
        n = cp / R

        T01 = T1 + Vm1**2 / (2.0 * cp)
        a0 = math.sqrt(gamma * R * T01)

        W1s = state.W1s
        T0R1s = T1 + W1s**2 / (2.0 * cp)
        P0R1s = P1 * (T0R1s / T1) ** n
        rho0R1s = P0R1s / (R * T0R1s)
        a0R1s = math.sqrt(gamma * R * T0R1s)

        A1 = state.mdot / (state.rho1 * Vm1)  # exact, from LE continuity
        # The REAL throat, when the caller supplies it (Impeller.throat_area). The old
        # A_th := A1 stand-in overstates the HECC inducer throat by 53% and makes this
        # loss structurally inert -- the model could not choke (debt D5).
        A_th = A1 if state.A_th is None else state.A_th
        beta_1xi = math.atan2(U1, Vm1)
        # Kovář et al. 2021 Eq. (21) (R7, docs/centrifugal/20-review-fixes.md): the sqrt was
        # dropped here previously -- restored. Reduces to sqrt(sin(beta_1xi)) iff A_th==A1.
        Cr = math.sqrt(A1 * math.sin(beta_1xi) / A_th)

        bracket = (2.0 + (gamma - 1.0) * (U1s / a0) ** 2) / (gamma + 1.0)
        A_star_th = (
            state.mdot / (rho0R1s * a0R1s) * bracket ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))
        )

        x = self.onset_scale * (
            self.onset_threshold - (Cr * A_th) / A_star_th
        )  # reciprocal ratio -- see docstring
        if x <= 0.0:
            return 0.0
        W1xi = math.hypot(Vm1, U1)
        # W1xi SQUARED. The transcription in research/notes 02 S2a (and the version this
        # file shipped until the throat was wired) reads `W1xi * (...)/2`, which is
        # DIMENSIONALLY IMPOSSIBLE: W is a velocity (m/s) and the bracket is
        # dimensionless, so it returns m/s, not J/kg = m^2/s^2. Every other loss in this
        # module is 1/2 * V^2. docs/PHYSICS-RULES.md rule 1.
        #
        # It survived precisely because it was STRUCTURALLY ZERO (A_th := A1 -> x < 0
        # always), so it never produced a number anyone could dimension-check -- the same
        # way the upstream defects survived. Wiring the real throat made it produce 367
        # J/kg where ~18,000 is needed at choke: 50x low, and dimensionally impossible.
        dh = W1xi**2 * (self.k_low * x + x**self.high_power) / 2.0
        return max(dh, 0.0)


def throat_relative_velocity(state: ImpellerLossState) -> float:
    """Inducer-throat relative velocity, ``Wth`` -- continuity + rothalpy + isentropic
    LE -> throat, the same system Kosuge, Ito & Nakanishi (1982), *ASME J. Eng. Power*
    104(4):782-787, pose as "three nonlinear algebraic equations".

    ⚠ **R8, docs/centrifugal/20-review-fixes.md -- PROVENANCE, not physics.** This
    closure had cited "Meroni, Zuhlsdorf, Elmegaard & Haglind (2018), *Applied
    Energy* 232:139-156, Eqs. (3)-(5)" for the three equations below, and Kosuge
    (1982) "via Meroni 2018 / Kovář et al. 2021 Eqs. (23)-(24)". **This project does
    not hold the Meroni 2018 paper** (``research/papers/`` has no Meroni PDF or
    text extraction, and no ``research/notes/`` file mentions it), and Kosuge's own
    primary is separately PAYWALLED and NOT HELD -- so a claim was being cited
    "via" a source that is *itself* unheld, and specific equation numbers (3)-(5)
    were asserted for a paper never opened. **The equation numbers below are
    therefore marked `# UNVERIFIED -- paper NOT HELD`, not removed**: the physics
    itself (continuity, rothalpy conservation across a rotor per
    docs/PHYSICS-RULES.md rule 3, isentropic LE->throat) is elementary and is
    independently verified by this module's OWN subsonic-branch bisection agreeing
    to 1e-6 with an independently-coded ``brentq`` solve of the identical system
    (``tests/unit/test_inducer_throat.py``) -- that agreement, not the Meroni
    citation, is what this closure's correctness rests on.

        Wth    = mdot / (rho_th * A_th)                    # continuity at the throat  # UNVERIFIED -- paper NOT HELD (cited as Meroni 2018 Eq. (3))
        h_th   = h01 - 0.5*Wth**2 + 0.5*U1_rms**2           # rothalpy, LE -> throat    # UNVERIFIED -- paper NOT HELD (cited as Meroni 2018 Eq. (4))
        s_th   = s01                                        # isentropic, LE -> throat  # UNVERIFIED -- paper NOT HELD (cited as Meroni 2018 Eq. (5))

    For a constant-cp ideal gas, with ``U_th`` taken at the LE RMS radius (so the
    rothalpy's ``-U**2/2`` term is UNCHANGED LE -> throat -- docs/PHYSICS-RULES.md
    rule 3 would otherwise require carrying a separate ``U_th``), this collapses to a
    1-D root solve for ``Wth``:

        T0R1   = T1 + W1xi**2/(2*cp)              # relative stagnation temp at the LE
        T_th   = T0R1 - Wth**2/(2*cp)
        rho_th = rho1 * (T_th/T1)**(1/(gamma-1))  # isentropic: s_th = s01 = s1
        residual(Wth) = rho_th*Wth*A_th - mdot

    Two roots exist (the throat mass flux rises with ``Wth`` up to ``Ma_th = 1``, then
    falls). This takes the SUBSONIC branch by bisection from the low side -- the SAME
    discipline as ``solve_vm_from_continuity``/``_solve_for_meridional_velocity``
    (``turbodesign.centrifugal.solver``/``diffusion``; docs/PHYSICS-RULES.md rule 5:
    bound the search on the local Mach number, never hand-roll a Newton step that can
    jump to the supersonic root).

    **RETRACTED, slice S3 (docs/centrifugal/15-experiments-prereg.md E2).** The plan
    pre-registered Wth = 217.4 +- 2 m/s against W1xi = 216.4 m/s -- "the throat is
    matched to 0.3%." The honest closure, using ``state.U1`` (the RMS-radius blade
    speed this module has used for "W1xi" since slice S1, not the arithmetic mean the
    pre-registered figure turns out to reconstruct from), gives Wth = 202.3 m/s
    against W1xi = 230.06 m/s at the HECC design point -- the throat is 7.9% LARGER
    than the LE relative-flow area, NOT matched (verified twice: this module's own
    bisection and an independently-coded ``brentq`` solve of the same system agree to
    1e-6; see ``tests/unit/test_inducer_throat.py``, which the shipped docstring used
    to contradict).

    FALLBACK: when ``state.A_th`` is ``None``, returns ``W1xi`` -- mirroring
    :class:`ImpellerChokeAungier`'s documented A_th-less fallback, so a caller that
    has not supplied a throat area gets an inert (zero) entrance-diffusion loss
    instead of a crash. This is NOT "the same approximation" the two classes used to
    share: :class:`ImpellerChokeAungier` was given a real throat and no longer falls
    back to this stand-in in practice; it survives here only for callers (bare unit
    fixtures) that never populate ``A_th`` at all.
    """
    W1xi = math.hypot(state.Vm1, state.U1)
    if state.A_th is None:
        return W1xi

    fluid = state.fluid
    gamma, cp, R = fluid.gamma, fluid.cp, fluid.R
    T1, rho1, A_th, mdot = state.T1, state.rho1, state.A_th, state.mdot

    T0R1 = T1 + W1xi**2 / (2.0 * cp)

    def flux(W: float) -> float:
        T_th = T0R1 - W * W / (2.0 * cp)
        if T_th <= 1.0:
            return -1.0
        rho_th = rho1 * (T_th / T1) ** (1.0 / (gamma - 1.0))
        return rho_th * W * A_th - mdot

    # Sonic (Ma_th = 1) relative velocity bounds the subsonic branch -- the same
    # closed form as _sonic_meridional_velocity(T0_eff, k=1, fluid) in solver.py, not
    # imported here (solver.py imports FROM this module, not the reverse).
    W_star = math.sqrt(2.0 * gamma * R * T0R1 / (gamma + 1.0))
    lo, hi = 1e-6, W_star
    if flux(lo) * flux(hi) > 0.0:
        raise ValueError(
            f"no subsonic inducer-throat solution for mdot={mdot:.4f} kg/s through "
            f"A_th={A_th:.6f} m2 (the throat is choked at this LE state, or A_th is "
            "too small)"
        )
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if flux(lo) * flux(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def kosuge_stall_ratio(state: ImpellerLossState) -> float:
    """Kosuge, Ito & Nakanishi (1982) inducer-stall criterion, ``W1s/Wth`` -- cited
    VIA Meroni 2018 / Kovář et al. 2021 Eqs. (23)-(24) (Kosuge's own primary is
    PAYWALLED and NOT HELD; the 0.5/1.75 constants below are therefore CITED-ONLY,
    not verified against the primary).

    Stall is predicted when this ratio exceeds 1.75 (a 0.5 "recovery factor" scales
    the corresponding loss in Kosuge's own model, which this module does not
    implement). **Diagnostic only** -- docs/centrifugal/17-tdd-plan.md slice S3: FOR
    FREE from the same throat solve as :func:`throat_relative_velocity`, NOT wired
    into this or any other loss.

    HECC, design point: W1s/Wth = 1.33 -- comfortably below the 1.75 stall
    threshold; the ratio is predicted to cross 1.75 at roughly 72% of design mass
    flow (docs/centrifugal/15-experiments-prereg.md E2).
    """
    return state.W1s / throat_relative_velocity(state)


def inducer_choke_mdot(state: ImpellerLossState) -> float:
    """Diagnostic: the mass flow at which the inducer throat first reaches
    ``Ma_th = 1`` -- i.e. where :func:`throat_relative_velocity`'s subsonic and
    supersonic branches merge -- holding the CURRENT LE relative-stagnation state
    (``T0R1``, ``T1``, ``rho1``) fixed. FOR FREE from that same solve
    (docs/centrifugal/17-tdd-plan.md slice S3); NOT wired into any loss.
    """
    if state.A_th is None:
        raise ValueError("inducer_choke_mdot needs a throat area (Impeller.throat_area)")
    fluid = state.fluid
    gamma, cp, R = fluid.gamma, fluid.cp, fluid.R
    W1xi = math.hypot(state.Vm1, state.U1)
    T0R1 = state.T1 + W1xi**2 / (2.0 * cp)
    W_star = math.sqrt(2.0 * gamma * R * T0R1 / (gamma + 1.0))
    T_star = T0R1 - W_star**2 / (2.0 * cp)
    rho_star = state.rho1 * (T_star / state.T1) ** (1.0 / (gamma - 1.0))
    return rho_star * W_star * state.A_th


@dataclass(frozen=True)
class ImpellerEntranceDiffusionAungier:
    """Aungier entrance-diffusion loss -- Kovář et al. 2021 Eqs. (22)-(24) (verified
    transcription of Aungier 2000). NOT in Oh 1997's canonical set; required here for
    the SAME reason as :class:`ImpellerChokeAungier` (Kovář et al. 2021 S7.1's documented
    near-choke over-prediction). INTERNAL: corrects incidence models under-predicting
    loss at positive incidence, where the dominant effect is diffusion from the LE to
    the inducer throat (Kovář et al. 2021 section 3.1.3 prose: "the incidence equations
    often underestimate the entrance loss at positive incidence angles ... the
    entrance diffusion loss ... corrects this").

        Δh_edf = max(0, 0.4*(W1xi - Wth)^2 - Δh_inc),   applied ONLY when W1xi > Wth

    **UN-PHANTOMED, slice S3 (docs/centrifugal/17-tdd-plan.md).** This class used to
    take ``Wth := W1xi`` -- "consistent with :class:`ImpellerChokeAungier`'s SAME
    approximation (throat == LE annulus)". That justification EXPIRED the moment
    ``ImpellerChokeAungier`` was given a real throat (``Impeller.throat_area =
    0.020525 m2``, NASA Appendix C) and started using it: this class was left behind,
    with ``state.A_th`` populated and sitting unused right next to it, making
    ``W1xi - Wth`` identically 0 and this loss identically 0 at EVERY operating point,
    forever -- the same silent-zero defect (docs/PHYSICS-RULES.md defect D) this
    project keeps finding. A term that cannot fire cannot be caught being wrong.

    ``Wth`` now comes from :func:`throat_relative_velocity` -- continuity + rothalpy
    + isentropic LE -> throat, needing only ``A_th`` (R8, docs/centrifugal/
    20-review-fixes.md: that function's docstring previously cited specific Meroni
    et al. (2018) equation numbers for a paper this project does not hold; see its
    docstring for the correction -- the physics is unchanged, the citation is not
    asserted past what is verified).

    **PRE-REGISTERED NULL, REFUTED (docs/centrifugal/15-experiments-prereg.md E2).**
    The plan predicted the design point would be UNCHANGED: HECC's inducer throat
    "matched" to its design-point relative flow to 0.3%, so
    ``0.4*(216.4-217.4)^2 = 0.4 J/kg``, minus ``Δh_inc`` (61 J/kg), clamped to zero.
    That prediction used a mean-radius ``U1`` this module has never used for "W1xi";
    with the RMS-radius ``state.U1`` every OTHER class in this module already reads
    (:class:`ImpellerIncidenceConrad`, :class:`ImpellerMixingAungier`,
    :class:`ImpellerChokeAungier`), the throat (``A_th = 0.020525 m2``) is genuinely
    7.9% LARGER than the LE relative-flow area (``A1*cos(beta1) = 0.01903 m2``), so
    the passage widens, not narrows, LE->throat: Wth = 202.3 m/s vs W1xi = 230.06
    m/s, and this loss fires **248.3 J/kg AT THE DESIGN POINT** -- comfortably over
    the plan's own 200 J/kg refutation threshold, not a null. It is NOT phantom any
    more, and it does not clamp to zero at design. It also fires OFF-DESIGN, growing
    monotonically as flow falls: 665.7 J/kg at 90% flow, 1155.2 at 80%, 1708.9 at 70%
    (Kosuge's stall branch active there -- ``W1s/Wth > 1.75`` -- see
    :func:`kosuge_stall_ratio`). See ``tests/unit/test_inducer_throat.py`` for the
    full accounting (including how the pre-registered figure reconstructs almost
    exactly under the mean-radius convention) -- the shipped docstring used to
    contradict that same test suite.

    **THE SIGN GUARD, added this slice, `# UNVERIFIED`.** Kovář et al. Eq. (22) as printed is
    unguarded on sign and therefore fires on ACCELERATION too (``Wth > W1xi``, i.e.
    toward choke) -- checked numerically, ~1700 J/kg of spurious "diffusion" loss
    above the design mdot, in a region :class:`ImpellerChokeAungier` already owns: a
    double count. This class applies the loss ONLY when ``W1xi > Wth`` (genuine
    LE->throat DECELERATION). Whether Aungier's own primary (Aungier 1995, *J.
    Turbomach.* 117(3):360-366 -- PAYWALLED, NOT HELD) restricts it so is
    UNVERIFIED; see data/coefficients.md.

    **``- Δh_inc`` stays -- a DOUBLE-COUNT GUARD, not an addition.** Both losses model
    the SAME LE->throat control volume: ``Δh_inc + Δh_edf = max(Δh_inc,
    0.4*(W1xi-Wth)^2)`` (Kovář et al. 2021 section 3.1.3 prose, quoted above). This class
    instantiates its OWN ``ImpellerIncidenceConrad()`` (default ``f_inc=0.6``) to
    compute that guard term -- FLAGGED, not fixed: if a caller ever configures a
    ``LossSet`` with an ``ImpellerIncidenceConrad`` whose ``f_inc`` has been
    overridden away from the default, this guard silently uses the WRONG (default)
    ``f_inc`` and desynchronizes from the incidence loss actually being reported. A
    one-line follow-up if it ever bites: thread the configured incidence model in
    instead of constructing a fresh default one.

    Coefficients: 0.4 (verified); 0.5 and 1.75 for Kosuge's stall criterion --
    CITED-ONLY diagnostic (:func:`kosuge_stall_ratio`), not wired into this or any
    other loss.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    k: float = FROZEN[
        "impeller.entrance_diffusion.k"
    ]  # Aungier (2000) via Kovář et al. 2021 Eqs. (22)-(24)

    def delta_h(self, state: ImpellerLossState) -> float:
        W1xi = math.hypot(state.Vm1, state.U1)
        Wth = throat_relative_velocity(state)
        if W1xi <= Wth:
            # Sign guard: only a genuine LE->throat DECELERATION diffuses -- firing on
            # acceleration (Wth > W1xi) would double-count against
            # ImpellerChokeAungier's near-choke physics. # UNVERIFIED against
            # Aungier's primary (paywalled) -- see class docstring.
            return 0.0
        # FLAG (not fixed): a fresh, default-f_inc incidence model -- may desync from
        # a caller-configured ImpellerIncidenceConrad. See class docstring.
        dh_inc = ImpellerIncidenceConrad().delta_h(state)
        dh = self.k * (W1xi - Wth) ** 2 - dh_inc
        return max(dh, 0.0)


# -------------------------------------------------------------------------- parasitic


@dataclass(frozen=True)
class ImpellerDiscFrictionDaily:
    """Daily & Nece (1960) disc-friction (windage) loss -- Oh, Yoon & Chung (1997)'s own
    choice (research notes S8a; Kovář et al. 2021 Eqs. (52)-(53), a verified transcription
    of Daily, J.W. & Nece, R.E., "Chamber Dimension Effects on Induced Flow and Frictional
    Resistance of Enclosed Rotating Disks", ASME J. Basic Eng. 82(1), 1960). PARASITIC:
    torque on the impeller BACK FACE (the unshrouded side facing the casing, outside the
    primary flow passage) that consumes shaft work and heats the fluid with no static-
    pressure benefit.

        Δh_df = 0.25 * Kf * rho_bar * D2^2 * U2^3 / (4 * mdot),   rho_bar = (rho1+rho2)/2

        Kf = 3.7 * G^0.1 / sqrt(Re2)     Re2 <= 3e4   (laminar)
           = 0.102 * G^0.1 / Re2^0.2     Re2 >  3e4   (turbulent)

        Re2 = U2 * D2 / nu2,   nu2 = mu / rho2,   G = 2*eps/D2 (back-face axial-gap ratio)

    **THE "TRANSCRIPTION ARTEFACT" CAVEAT ON THE LEADING 0.25 IS WITHDRAWN.** Notes S8a --
    and, until now, this docstring and ``data/coefficients.md`` -- flagged the leading
    ``0.25 * (.../4)`` as a "possible transcription artefact in Kovář et al., reproduced
    exactly". **It is not an artefact. It is a unit conversion, and the code is right.**
    Oh, Yoon & Chung (1997), Table 6, p. 336 prints the disc-friction row as

        Δh_df = f_df * rho_bar * r2**2 * U2**3 / (4 * mdot)     <- r2 SQUARED, not D2

    and since ``D2**2 = 4*r2**2``, writing it over D2 costs exactly a factor 1/4: the 0.25
    IS the D2 -> r2 conversion, algebraically identical to Oh. Checked numerically against
    Oh's own form: this class returns **1.058x** Oh's value, and the whole of that residual
    is Kovář's clearance term inside Kf, not the prefactor. **The VERBATIM stamp in
    data/coefficients.md is CORRECT and stays.** The laminar/turbulent
    transition Reynolds number is ALSO inconsistent across the literature (Kovář et al.: 3e4;
    radcomp/CIMdes: 3e5) -- flagged, not resolved; HECC's Re2 is ~3e7 either way (see
    below), so the choice is immaterial here.

    ``G`` (the back-face axial-clearance ratio) is genuinely unmeasured for this
    machine: neither ``Impeller`` nor ``MeridionalPath`` carry a back-face cavity model.
    Daily & Nece's own experiments spanned ``G = 0.0127-0.217``; the default below
    (0.05) is the order-of-magnitude MIDPOINT of that TESTED range, not a value read off
    HECC's drawings -- **# UNVERIFIED for this specific machine.** It is not a hidden
    tuning knob, though: the exponent on G is 0.1, so at the HECC design point psi moves
    only 0.7803 -> 0.7819 (0.002, about 5% of the Slice 4 psi tolerance) across the
    ENTIRE Daily & Nece tested range -- checked numerically, not asserted by a test.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "parasitic"
    clearance_ratio: float = FROZEN[
        "impeller.disc_friction.back_face_clearance_ratio"
    ]  # G = 2*eps/D2; UNVERIFIED for HECC -- see class docstring
    re_transition: float = FROZEN[
        "impeller.disc_friction.re_transition"
    ]  # laminar/turbulent Re2 boundary
    laminar_c: float = FROZEN[
        "impeller.disc_friction.laminar_c"
    ]  # Daily & Nece (1960) via Kovář et al. 2021 Eq. (52)
    turbulent_c: float = FROZEN[
        "impeller.disc_friction.turbulent_c"
    ]  # Daily & Nece (1960) via Kovář et al. 2021 Eq. (53)
    g_exponent: float = FROZEN["impeller.disc_friction.g_exponent"]
    re_exponent: float = FROZEN["impeller.disc_friction.re_exponent"]
    leading_factor: float = FROZEN["impeller.disc_friction.leading_factor"]
    denominator: float = FROZEN["impeller.disc_friction.denominator"]
    mu: Optional[float] = None  # Pa*s; None -> state.fluid.mu (ISA constant, state.py)

    def delta_h(self, state: ImpellerLossState) -> float:
        assert state.rho2 is not None, (
            "ImpellerDiscFrictionDaily needs the impeller-exit density, which is only "
            "known once the internal-loss/continuity solve has converged -- evaluate "
            "parasitic losses AFTER internal losses, not inside the same trial step "
            "(see turbodesign.centrifugal.losses module docstring)."
        )
        D2 = 2.0 * state.r2
        U2 = state.U2
        mu = self.mu if self.mu is not None else state.fluid.mu
        rho_bar = 0.5 * (state.rho1 + state.rho2)

        Re2 = max(abs(U2) * D2 * state.rho2 / mu, 1.0)
        G = self.clearance_ratio
        if Re2 <= self.re_transition:
            Kf = self.laminar_c * G**self.g_exponent / math.sqrt(Re2)
        else:
            Kf = self.turbulent_c * G**self.g_exponent / Re2**self.re_exponent

        dh = self.leading_factor * Kf * rho_bar * D2**2 * U2**3 / (self.denominator * state.mdot)
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerRecirculationOh:
    """Oh, Yoon & Chung (1997) recirculation loss -- "the only genuinely Oh-original
    equation in the Oh set" (research notes 02 S9b; Kovář et al. 2021 Eq. (49),
    verified; confirmed independently by CIMdes ``loss.py::RecirculationLoss``).
    PARASITIC: backflow from the diffuser into the impeller tip, growing sharply with
    exit swirl -- this is what shapes the surge-side efficiency rollover.

        Δh_rc = 8e-5 * sinh(3.5*alpha2^3) * D_f^2 * U2^2,     alpha2 in RADIANS

    Coefficients: 8e-5; 3.5; the cubic power inside sinh (verified). ``D_f`` from
    :func:`_diffusion_factor` (the same one the blade-loading loss uses -- research
    notes 01 S4.2: "D_f ... drives #4 and #9").

    ⚠ **Known bug, avoided.** Reading ``turbodesign/loss/compressor/otac.py`` line
    669 (``ImpellerRecirculationOh``, upstream's translation of this same
    correlation) turned up ``np.sinh(3.5 * (np.radians(row.alpha2) ** 3))`` where
    ``row.alpha2`` is ALREADY in radians -- a SECOND ``np.radians()`` call shrinks a
    ~0.5 rad angle to ~0.009 rad, cubes it to ~7e-7, and the loss comes out about
    five orders of magnitude too small, silently (effectively zeroing it). This
    implementation computes ``alpha2`` directly from ``atan2`` (radians by
    construction, the Python ``math`` module's convention) and applies NO further
    unit conversion.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "parasitic"
    k: float = FROZEN[
        "impeller.recirculation.k"
    ]  # Oh, Yoon & Chung (1997) via Kovář et al. 2021 Eq. (49)
    sinh_arg_coeff: float = FROZEN["impeller.recirculation.sinh_arg_coeff"]
    alpha2_power: float = FROZEN["impeller.recirculation.alpha2_power"]

    # Warn above this exit swirl angle. sinh(3.5*alpha2**3) is an EXPONENTIAL amplifier
    # and this term is the largest parasitic contributor by far (89.5% of parasitic work
    # at the HECC design point). Its argument grows as the CUBE of alpha2:
    #
    #     alpha2      3.5*alpha2**3      sinh(.)     vs the HECC design point
    #      60 deg          4.02             27.7          0.014x
    #      70 deg          6.39            296.5          0.153x
    #      76.3 deg        8.26           1931.3          1.000x   <- HECC design point
    #      80 deg          9.19           4877.3          2.525x
    #
    # so a 1-DEGREE error in alpha2 near the HECC point swings the dominant parasitic
    # term by 30-40% -- roughly the entire slice-5 psi gate. And alpha2 is not an input:
    # it is a COMPUTED output of the exit triangle (slip x backsweep x blockage x Cm2),
    # none of which is independently validated. docs/PHYSICS-RULES.md: "A correlation outside its
    # fitted range is an error, not an extrapolation. Guard and warn."
    #
    # 70 deg: # UNVERIFIED. Oh, Yoon & Chung (1997) state no validity range for alpha2,
    # and the paper's own impeller set is not tabulated in any source we hold. This
    # threshold is ENGINEERING JUDGEMENT, chosen as the angle beyond which sinh has
    # clearly left its near-linear region (sinh(x)/x = 46 at 70 deg, 234 at 76.3 deg) --
    # NOT a published bound. It is a tripwire, not a physical limit. Do not cite it as Oh's.
    alpha2_warn_deg: float = FROZEN["impeller.recirculation.alpha2_warn_deg"]

    def delta_h(self, state: ImpellerLossState) -> float:
        Df = _diffusion_factor(state)
        alpha2 = math.atan2(state.Vt2, state.Cm2)  # radians already -- no further
        # math.radians()/np.radians() call (see the otac.py bug noted above).

        if math.degrees(alpha2) > self.alpha2_warn_deg:
            x = self.sinh_arg_coeff * alpha2**self.alpha2_power
            # d(sinh x)/d(alpha2) / sinh x = 3*x/alpha2 * coth(x) ~ 3*x/alpha2 for x >> 1
            per_deg = 100.0 * (self.alpha2_power * x / alpha2) * math.radians(1.0) / math.tanh(x)
            warnings.warn(
                f"ImpellerRecirculationOh evaluated at alpha2 = {math.degrees(alpha2):.1f} deg "
                f"(> {self.alpha2_warn_deg} deg): sinh argument = {x:.2f}, sinh = "
                f"{math.sinh(x):.0f}. This is deep in the correlation's EXPONENTIAL "
                f"region, where a 1-degree error in alpha2 changes this loss by "
                f"~{per_deg:.0f}%. alpha2 is a computed output, not a measurement, and "
                f"this term dominates the parasitic budget -- treat the resulting work "
                f"factor as weakly constrained. (Threshold is UNVERIFIED engineering "
                f"judgement, not a published bound; Oh 1997 states no validity range.)",
                RuntimeWarning,
                stacklevel=2,
            )

        dh = (
            self.k
            * math.sinh(self.sinh_arg_coeff * alpha2**self.alpha2_power)
            * Df**2
            * state.U2**2
        )
        return max(dh, 0.0)


@dataclass(frozen=True)
class ImpellerLeakageAungier:
    """Aungier tip-leakage loss -- Oh, Yoon & Chung (1997)'s selection for leakage
    (research notes 02 S10a). PARASITIC: extra shaft work spent pumping the
    recirculated tip-leakage mass flow -- distinct from the clearance loss's (S6a)
    aerodynamic tip-vortex penalty.

        Δh_lk = mdot_cl * Ucl * U2 / (2*mdot)     # REPAIRED, E-R11 (C2)
        Ucl = 0.816 * sqrt(2*dP_cl / rho2)
        mdot_cl = (rho1+rho2) * Z * eps * Lm * Ucl / 2
        dP_cl = mdot*(D2*Cu2 - D1xi*Cu1xi) / (Z*Dbar*bbar*Lm)
        Dbar = (D1xi+D2)/2,   bbar = (b1+b2)/2,   Cu1xi = 0 (no inlet swirl)

    Coefficient: 0.816 (= sqrt(2/3), the gap discharge coefficient; verified,
    research notes 02 S6d). ``eps`` is ``Impeller.tip_clearance``.

    ========================= REPAIRED, E-R11 (C2). DEBT D6 CLOSED. =========================
    This term shipped for months as ``mdot_cl*Ucl/(2*mdot)`` -- dimensionally **m/s**, not
    J/kg, one factor short of a specific energy. It was left BROKEN, deliberately, under
    one instruction: *"DO NOT GUESS THE MISSING FACTOR."* That instruction was RIGHT, and
    it is now VINDICATED BY THE PRIMARIES.

    *** THE MISSING FACTOR IS U2. IT IS NOT Ucl. ***

    TWO INDEPENDENT PRIMARIES, BOTH NOW HELD, AGREE:

      * Aungier (1995) Eq. (10), OPEN impellers:
            I_DF + I_L = C_MD*rho2*U2*r2**2/(2*mdot)  +  mdot_cl*Ucl/(2*mdot*U2)
        and ``I`` is a WORK-INPUT COEFFICIENT, I == dh/U2**2 (his Eq. 1 + Nomenclature),
        so   dh_L = I_L * U2**2 = mdot_cl * Ucl * U2 / (2*mdot).
      * Oh, Yoon & Chung (1997) Table 6 prints, independently:
            dh_lk = mdot_cl * Ucl * U2 / (2*mdot)

    THE OLD DOCSTRING CALLED ``* Ucl`` "the structurally obvious candidate" FOR THE
    MISSING FACTOR. IT WAS THE WRONG ONE. Both ``*Ucl`` and ``*U2`` restore J/kg. Only
    ``*U2`` is the physics -- and they are not close: U2 = 492.4 m/s on HECC, while Ucl is
    a gap throughflow velocity an order of magnitude smaller.

        DIMENSIONAL ANALYSIS DOES NOT DETERMINE A MISSING FACTOR.
        TWO DIFFERENT VELOCITIES RESTORE THE UNITS; ONLY ONE IS THE PHYSICS.

    AND THE SOURCE WAS ON OUR SHELF THE WHOLE TIME. This docstring previously said the
    citation was "FABRICATED / NOT FOUND" and that *"Aungier (2000)'s own leakage section
    MUST BE OBTAINED before this formula is repaired."* IT DID NOT NEED TO BE. The correct
    formula is in **Table 6 of Oh, Yoon & Chung (1997)** -- the paper this loss set is
    NAMED AFTER, held in research/papers/ from the start, and already cited as the
    authority for EIGHT OTHER TERMS IN THIS FILE. We searched for an Aungier term in an
    Aungier paper and never opened the Oh paper that transcribes it.

        A SHELF IS NOT A SEARCH. HOLDING A SOURCE IS NOT READING IT.

    (The Kovář et al. 2021 S3.1.6 attribution remains WRONG and is not restored: that
    section's Eqs. (38)-(41) are the *clearance* loss, a different term.)

    PARASITIC (rule 6): it ADDS WORK, it does NOT destroy pressure. Measured on adoption:
    dh_lk = 936 J/kg on HECC (was ~1.91, i.e. really 1.91 m/s), and stage PR moved by
    0.020 pp -- rule 6 holds, and the solver does not leak parasitic work into pressure.
    See docs/centrifugal/40-aungier-primary-result.md S3 and 43-oh-primary-result.md S4.
    =========================================================================================

    ``Lm`` (impeller MERIDIONAL length, LE to TE along the flow path -- A DIFFERENT
    QUANTITY from the through-blade camberline length ``L_tilde`` the skin-friction and
    mixing losses use) is ``state.L_meridional`` (``Impeller.l_main_m``): NASA's own
    measured meridional camberline length for HECC (~0.2099 m, already sourced and used
    by ``Impeller.Z_eff``), used directly when supplied.

    **CLOSED, slice S1 (docs/centrifugal/17-tdd-plan.md).** Before this slice, ``Lm``
    was proxied by ``_camberline_length(state)`` -- the THROUGH-BLADE length -- "an
    unsourced closure doing double duty" (docs/centrifugal/
    12-model-as-implemented.md S7.2): the same closure fed skin friction, mixing, AND
    leakage, three consumers that structurally want two different quantities (a
    through-blade length for the first two, a meridional length for this one). Absent
    a measured ``L_meridional``, this falls back to ``_camberline_length(state)`` --
    the SAME proxy as before, for a caller that has not supplied either measured
    length (matching the ``Z_le``/``A_th`` "old behaviour by default" pattern
    elsewhere in :class:`ImpellerLossState`).

    ``b1`` (inlet blade height) is taken as ``r1s - r1h`` -- the LE annulus span, an
    EXACT geometric quantity (not an approximation of a different one).
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "parasitic"
    discharge_coeff: float = FROZEN[
        "impeller.leakage.discharge_coeff"
    ]  # Aungier gap discharge coeff, sqrt(2/3)

    def delta_h(self, state: ImpellerLossState) -> float:
        assert state.rho2 is not None, (
            "ImpellerLeakageAungier needs the converged impeller-exit density -- "
            "evaluate parasitic losses AFTER internal losses/continuity have "
            "converged, not inside the same trial step."
        )
        r1h, r1s, r2, b2 = state.r1h, state.r1s, state.r2, state.b2
        Z = state.Z
        Cu2 = state.Vt2
        D2 = 2.0 * r2
        r1xi = math.sqrt(0.5 * (r1h**2 + r1s**2))
        D1xi = 2.0 * r1xi
        b1 = r1s - r1h
        # MERIDIONAL length, NOT the through-blade camberline length -- a genuinely
        # different quantity (see class docstring). Falls back to the old
        # (through-blade) proxy only when the caller has supplied neither.
        Lm = state.L_meridional if state.L_meridional is not None else _camberline_length(state)

        Dbar = (D1xi + D2) / 2.0
        bbar = (b1 + b2) / 2.0
        dP_cl = state.mdot * (D2 * Cu2) / (Z * Dbar * bbar * Lm)  # Cu1xi = 0

        Ucl = self.discharge_coeff * math.sqrt(max(2.0 * dP_cl / state.rho2, 0.0))
        mdot_cl = (state.rho1 + state.rho2) * Z * state.tip_clearance * Lm * Ucl / 2.0

        # REPAIRED, E-R11 (C2). The missing factor is U2 -- NOT Ucl.
        #
        #   Aungier (1995) Eq. (10), OPEN impellers:  I_L = mdot_cl*Ucl/(2*mdot*U2),
        #   and I == dh/U2**2, so  dh_L = mdot_cl * Ucl * U2 / (2*mdot)   [J/kg].
        #   Oh, Yoon & Chung (1997) Table 6 prints the SAME formula, independently.
        #   TWO PRIMARIES AGREE.
        #
        # This class previously shipped `mdot_cl*Ucl/(2*mdot)` -- dimensionally m/s --
        # and its docstring named `* Ucl` "the structurally obvious candidate" for the
        # missing factor. IT WAS THE WRONG ONE. Both restore J/kg; only U2 is the
        # physics, and on HECC they differ by ~an order of magnitude (U2 = 492.4 m/s).
        # DIMENSIONAL ANALYSIS DOES NOT DETERMINE A MISSING FACTOR. Refusing to guess
        # was right, and the primaries settled it -- see docs/centrifugal/40- and 43-.
        dh = mdot_cl * Ucl * state.U2 / (2.0 * state.mdot)
        return max(dh, 0.0)


# ----------------------------------------------------------------------- the Oh 1997 set


@dataclass(frozen=True)
class OhLossSet(LossSet):
    """The Oh, Yoon & Chung (1997) canonical loss set, plus Aungier's choke and
    entrance-diffusion additions (research notes 02 S16-S17; Kovář et al. 2021 Table
    3, "Set Oh"). See ``tests/unit/test_slice5_full_loss_set.py``'s module docstring
    for why the two Aungier additions are not optional here (Oh's set alone "largely
    overestimat[es] efficiency near choke point", Kovář et al. 2021 S7.1).

    INTERNAL (destroy P02, add no work): incidence (Conrad), blade loading (Coppage),
    skin friction (Jansen), clearance (Jansen), mixing (Aungier -- see
    :class:`ImpellerMixingAungier`'s docstring for why this deviates from Oh's
    literal Johnston & Dean choice), choke (Aungier), entrance diffusion (Aungier).

    PARASITIC (add work, add no pressure): disc friction (Daily & Nece), recirculation
    (Oh), leakage (Aungier).

    Every coefficient inside every model below is frozen at its cited published
    value -- this class takes no constructor arguments of its own (the individual
    classes still keep their own calibration knobs, e.g.
    ``ImpellerSkinFrictionJansen(Cf=...)``, for callers who explicitly want to vary
    one). ``OhLossSet()`` exists to be the ONE, un-tunable "reference Oh set"
    configuration this project validates psi/eta against.

    ⚠️ **THIS SET IS A HYBRID THAT NO PUBLISHED PAPER EVER VALIDATED AS A SET.** Oh's own
    nine terms (E2, docs/centrifugal/43-oh-primary-result.md) contain **NO choke loss and
    NO entrance-diffusion loss**, and his mixing loss is **Johnston & Dean's, not
    Aungier's**. We swapped one term and added two. Neither the substitution nor the union
    was validated by Oh, by Aungier, or by anyone else. **Say so, every time.**

    ``apply_head_loss_correction`` defaults to **True** here: Aungier's ``f_c`` (Eq. 32)
    was ADOPTED in E-R11 with PI sign-off. Pass ``False`` to reproduce the PRE-ADOPTION
    frozen set exactly -- that configuration remains the published CONTROL.
    """

    models: Sequence[LossModel] = field(
        default_factory=lambda: (
            ImpellerIncidenceConrad(),
            ImpellerBladeLoadingCoppage(),
            ImpellerSkinFrictionJansen(),
            ImpellerClearanceJansen(),
            ImpellerMixingAungier(),
            ImpellerChokeAungier(),
            ImpellerEntranceDiffusionAungier(),
            ImpellerDiscFrictionDaily(),
            ImpellerRecirculationOh(),
            ImpellerLeakageAungier(),
        )
    )
    apply_head_loss_correction: bool = True


__all__ = [
    "EvaluatedLosses",
    "ImpellerBladeLoadingCoppage",
    "ImpellerChokeAungier",
    "ImpellerClearanceJansen",
    "ImpellerDiscFrictionDaily",
    "ImpellerEntranceDiffusionAungier",
    "ImpellerIncidenceConrad",
    "ImpellerLeakageAungier",
    "ImpellerLossState",
    "ImpellerMixingAungier",
    "ImpellerRecirculationOh",
    "ImpellerSkinFrictionJansen",
    "LossModel",
    "LossSet",
    "NO_LOSSES",
    "OhLossSet",
]
