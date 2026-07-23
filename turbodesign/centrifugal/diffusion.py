"""The diffusion system: vaneless space, vaned diffuser, exit guide vanes.

WHY THIS EXISTS (docs/centrifugal/07-diffusion-system-design.md)
---------------------------------------------------------------
Without it the model stops at the impeller trailing edge, so its PR (5.40) and eta (0.90)
are impeller-only and comparable to NOTHING -- NASA quotes STAGE values (4.6847, 0.8553).
That left exactly ONE constrained output, psi, and an overfitting audit showed what that
permits: psi is a one-parameter family in backsweep, and along a continuous
(backsweep, Z_eff) ridge psi = 0.81 is hit exactly while PR moves only 5.20-5.22 and eta
only 0.909-0.912. No model output could distinguish those parameter sets. The fit was
formally UNIDENTIFIABLE.

A diffuser is therefore not "one more component": it is what makes the model FALSIFIABLE.
With a stage exit, PR, eta and psi become checkable SIMULTANEOUSLY, and a one-parameter
family cannot satisfy three independent constraints.

FRAME DISCIPLINE (docs/PHYSICS-RULES.md rule 2)
-----------------------------------
These are STATIONARY rows: U = 0, so the relative frame IS the absolute frame and rothalpy
reduces to enthalpy. Consequently:

    T0 IS CONSERVED ACROSS EVERY COMPONENT HERE. No work is done.

That is a strong, cheap invariant, and every component in this module is tested against
it. A stationary component that changes T0 has a bug -- there is no physics that permits
it. Losses show up in P0 only.

SLIP MODELS ARE FORBIDDEN HERE (docs/PHYSICS-RULES.md rule 7)
-------------------------------------------------
A slip model must NEVER be applied to a vaned diffuser. Slip is the centrifugal analogue
of deviation and its derivation rests on the Coriolis term; with omega = 0 that term
vanishes and the derivation collapses. A vaned diffuser is a CASCADE, not a rotor. Use
Carter's rule. :class:`VanedDiffuser` and :class:`ExitGuideVanes` accept no slip model and
there is no code path by which one could reach them.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Optional, Protocol, runtime_checkable

from .coefficients import FROZEN
from .state import Air, entropy


@dataclass(frozen=True)
class FlowState:
    """The absolute-frame flow state at a stationary station.

    Deliberately NOT ``ThermoState``: that carries U, W, P0R, T0R and rothalpy, none of
    which mean anything where U = 0. Carrying them would invite exactly the frame confusion
    docs/PHYSICS-RULES.md rule 2 warns about. A stationary station has an absolute state and nothing else.
    """

    P0: float  # absolute total pressure, Pa
    T0: float  # absolute total temperature, K -- CONSERVED across every component here
    P: float  # static pressure, Pa
    T: float  # static temperature, K
    rho: float  # static density, kg/m3
    Vm: float  # meridional velocity, m/s
    Vt: float  # tangential velocity, m/s
    r: float  # radius, m
    b: float  # passage width, m
    s: float  # specific entropy relative to the stage inlet, J/(kg K)

    @property
    def V(self) -> float:
        return math.hypot(self.Vm, self.Vt)

    @property
    def alpha(self) -> float:
        """Absolute flow angle FROM MERIDIONAL, radians. 0 = purely meridional (no swirl),
        pi/2 = purely tangential. At a centrifugal impeller exit this is ~76 deg."""
        return math.atan2(self.Vt, self.Vm)

    @property
    def alpha_deg(self) -> float:
        return math.degrees(self.alpha)

    def mach(self, fluid: Air) -> float:
        return self.V / fluid.speed_of_sound(self.T)

    def mach_meridional(self, fluid: Air) -> float:
        """docs/PHYSICS-RULES.md rule 5: choke is on the MERIDIONAL Mach number. The ABSOLUTE Mach
        routinely exceeds 1 at a centrifugal impeller exit WITHOUT choking."""
        return self.Vm / fluid.speed_of_sound(self.T)


def state_from_totals(
    P0: float, T0: float, Vm: float, Vt: float, r: float, b: float, fluid: Air, s: float
) -> FlowState:
    """Build a static state from totals and the velocity triangle (isentropic statics)."""
    V2 = Vm * Vm + Vt * Vt
    T = T0 - V2 / (2.0 * fluid.cp)
    if T <= 0.0:
        raise ValueError(
            f"velocity {math.sqrt(V2):.1f} m/s exceeds the total enthalpy at T0={T0:.1f} K"
        )
    P = P0 * (T / T0) ** (fluid.gamma / (fluid.gamma - 1.0))
    return FlowState(P0=P0, T0=T0, P=P, T=T, rho=P / (fluid.R * T), Vm=Vm, Vt=Vt, r=r, b=b, s=s)


def solve_vm_from_continuity(
    P0: float, T0: float, Vt: float, area: float, mdot: float, fluid: Air
) -> float:
    """Meridional velocity satisfying rho(Vm)*Vm*area = mdot, at fixed P0, T0, Vt.

    Takes the SUBSONIC (meridional) root -- bisection from the low side. The mass flux
    rho*Vm rises with Vm up to Ma_m = 1 and falls after, so there are two roots; the
    subsonic one is the physical branch for an unchoked passage (docs/PHYSICS-RULES.md rule 5: the
    criterion is on the MERIDIONAL Mach number, not the absolute one, which is routinely
    supersonic here without choking).
    """

    def flux(Vm: float) -> float:
        V2 = Vm * Vm + Vt * Vt
        T = T0 - V2 / (2.0 * fluid.cp)
        if T <= 1.0:
            return -1.0
        P = P0 * (T / T0) ** (fluid.gamma / (fluid.gamma - 1.0))
        return (P / (fluid.R * T)) * Vm * area - mdot

    # Sonic meridional velocity bounds the subsonic branch.
    lo, hi = 1e-6, 1.0
    for _ in range(200):
        V2 = hi * hi + Vt * Vt
        T = T0 - V2 / (2.0 * fluid.cp)
        if T <= 1.0 or hi >= fluid.speed_of_sound(max(T, 1.0)):
            break
        hi *= 1.05
    if flux(lo) * flux(hi) > 0.0:
        raise ValueError(
            f"no subsonic solution for mdot={mdot:.4f} kg/s through area={area:.6f} m2 "
            f"(the passage is CHOKED at this state, or the area is too small)"
        )
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if flux(lo) * flux(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


@runtime_checkable
class StationaryComponent(Protocol):
    """A stationary flow-path component. Takes a state in, returns the state out.

    Every implementation MUST conserve T0 (U = 0 -> no work). Losses appear in P0 only.
    """

    name: str

    def solve(self, inlet: FlowState, mdot: float, fluid: Air) -> FlowState: ...


# --------------------------------------------------------------------- vaneless space


@dataclass(frozen=True)
class VanelessSpace:
    """The vaneless gap between the impeller TE and the diffuser vane LE (r2 -> r3).

    Short on HECC (radius ratio 1.072) but NOT a free pass: it is exactly where
    impeller/diffuser matching happens, and NASA blames HECC's efficiency shortfall partly
    on "an impeller/diffuser corrected flow-rate mismatch (and associated incidence
    levels)" (NASA/CR-2014-218114/REV1).

    PHYSICS. In the absence of friction, angular momentum is conserved exactly:

        r * Vt = const                      (free vortex)

    Friction is what actually sets the diffuser inlet flow angle -- and the inlet angle is
    what the vaned diffuser's incidence keys off, so it is first-order for the NEXT
    component even though its own pressure loss is modest. Integrating in r:

        d(r*Vt)/dr = -cf * r * V * Vt / (b * Vm)          # wall friction torque
        dP0/dr     = -cf * rho * V**3 / (2 * b * Vm)      # friction dissipation

    with continuity rho*Vm*2*pi*r*b*(1-blockage) = mdot closing Vm at each step, and
    T0 = const throughout (stationary row -- no work).

    Cf = 0.005: a standard smooth-wall vaneless-diffuser skin-friction coefficient
    (Japikse, "Centrifugal Compressor Design and Performance", 1996, Ch. 4; the same
    order as Aungier's 0.005 for a vaneless space). # UNVERIFIED against a primary
    tabulation -- recorded in data/coefficients.md. Cf = 0.0 must reduce this component
    EXACTLY to r*Vt = const and P0 = const, and there is a test asserting that.

    A PINCH IS SUPPORTED: b3 may differ from the impeller exit b2 (HECC pinches from
    0.609 in to 0.559 in, -8.2%). Passing b3 != b2 linearly tapers b across the gap.
    """

    r3: float  # vane leading-edge radius, m
    b3: float  # channel height at r3, m (may differ from b2 -- HECC has a pinch)
    blockage: float = 0.0
    Cf: float = FROZEN["diffuser.vaneless.cf"]  # see docstring; # UNVERIFIED
    n_steps: int = 50
    name: str = "VanelessSpace"

    def solve(self, inlet: FlowState, mdot: float, fluid: Air) -> FlowState:
        if self.r3 <= inlet.r:
            raise ValueError(f"vaneless space must diffuse outward: r3={self.r3} <= r2={inlet.r}")
        r, b2 = inlet.r, inlet.b
        rVt = inlet.r * inlet.Vt
        P0, T0 = inlet.P0, inlet.T0
        dr = (self.r3 - inlet.r) / self.n_steps

        for i in range(self.n_steps):
            b = b2 + (self.b3 - b2) * (i / self.n_steps)  # linear taper across the pinch
            area = 2.0 * math.pi * r * b * (1.0 - self.blockage)
            Vt = rVt / r
            Vm = solve_vm_from_continuity(P0, T0, Vt, area, mdot, fluid)
            V = math.hypot(Vm, Vt)
            T = T0 - V * V / (2.0 * fluid.cp)
            P = P0 * (T / T0) ** (fluid.gamma / (fluid.gamma - 1.0))
            rho = P / (fluid.R * T)

            # E-R3/D2 (docs/centrifugal/26-prereg-r3.md). These two lines MUST descend from
            # the same wall shear. With tau = Cf*rho*V^2/2 acting on BOTH endwalls, the
            # control-volume balance over an annulus of radial extent dr gives
            #
            #     d(r*Vt)/dr = -Cf * V * Vt * r / (Vm * b)      <- the torque line
            #     dP0/dr     = -Cf * rho * V**3   / (Vm * b)    <- the dissipation line
            #
            # The dissipation line previously carried an extra 1/2. That is not a Cf
            # convention: the RATIO of the two equations is rho*V^2/(r*Vt) and is INDEPENDENT
            # of how Cf is defined, so no definition of Cf can make both lines right. The
            # shipped ratio was exactly half of it -- the dissipation line applied HALF the
            # friction the swirl-decay line applied. Raised to the torque line's convention,
            # because the torque line is the one that matches the derivation.
            d_rVt = -self.Cf * r * V * Vt / (b * Vm) if self.Cf else 0.0
            dP0 = -self.Cf * rho * V**3 / (b * Vm) if self.Cf else 0.0
            rVt += d_rVt * dr
            P0 += dP0 * dr
            r += dr

        area3 = 2.0 * math.pi * self.r3 * self.b3 * (1.0 - self.blockage)
        Vt3 = rVt / self.r3
        Vm3 = solve_vm_from_continuity(P0, T0, Vt3, area3, mdot, fluid)
        out = state_from_totals(P0, T0, Vm3, Vt3, self.r3, self.b3, fluid, s=0.0)
        return _with_entropy(out, inlet, fluid)


# --------------------------------------------------------------------- cascade helpers


def carters_deviation(camber_deg: float, solidity: float) -> float:
    """Carter's rule: delta = m * theta / sqrt(solidity), degrees.

    ``theta`` is the blade camber angle (LE metal angle - TE metal angle) and ``solidity``
    is chord/pitch. ``m = 0.26`` is the standard circular-arc compressor-cascade value at
    moderate stagger (Carter, A.D.S. (1950), "The Low Speed Performance of Related
    Aerofoils in Cascade", ARC CP 29; as tabulated in Dixon & Hall, "Fluid Mechanics and
    Thermodynamics of Turbomachinery", 7th ed., Eq. 3.31 and Fig. 3.16).

    ⚠️ THIS -- NOT A SLIP MODEL -- IS THE DEVIATION RULE FOR A STATIONARY CASCADE
    (docs/PHYSICS-RULES.md rule 7). Slip rests on the Coriolis term; at omega = 0 that term vanishes
    and the derivation collapses. Applying a slip factor to a vaned diffuser is a physics
    error that returns a plausible number.

    m = 0.26: # UNVERIFIED at this precision -- Carter's m is a function of stagger and is
    read from a chart; 0.26 is the widely-quoted circular-arc value. See
    data/coefficients.md.
    """
    if solidity <= 0.0:
        raise ValueError(f"solidity must be > 0, got {solidity}")
    return FROZEN["diffuser.carter.m"] * camber_deg / math.sqrt(solidity)


def _with_entropy(out: FlowState, inlet: FlowState, fluid: Air) -> FlowState:
    """Carry entropy forward: ds = cp*ln(T/Tref) - R*ln(P/Pref), evaluated on the statics."""
    ds = entropy(out.T, out.P, inlet.T, inlet.P, fluid)
    return FlowState(
        P0=out.P0,
        T0=out.T0,
        P=out.P,
        T=out.T,
        rho=out.rho,
        Vm=out.Vm,
        Vt=out.Vt,
        r=out.r,
        b=out.b,
        s=inlet.s + ds,
    )


@dataclass(frozen=True)
class _Cascade:
    """Shared machinery for a stationary bladed row: incidence, deviation, friction.

    Both the vaned diffuser and the EGV are cascades. What differs is their geometry and
    which direction they turn the flow -- not the physics.
    """

    beta_le_deg: float  # vane LE metal angle, deg from meridional
    beta_te_deg: float  # vane TE metal angle, deg from meridional
    # FLOAT, not int: VanedDiffuser feeds this TWO different counts depending on the
    # station -- an integer TE/wetted-passage count (deviation, friction) and a
    # length-weighted fractional count (Z_eff_loading, diffusion factor / loading loss;
    # docs/centrifugal/17-tdd-plan.md slice S5), the same way Impeller.Z_eff is
    # fractional.
    n_vanes: float
    chord: float  # m
    Cf: float = FROZEN["diffuser.cascade.cf"]
    f_inc: float = FROZEN[
        "diffuser.cascade.f_inc"
    ]  # incidence loss coefficient, same basis as the impeller's Conrad
    blade_loading_k: float = FROZEN["diffuser.cascade.blade_loading_k"]
    df_warn_threshold: float = FROZEN["diffuser.cascade.df_warn_threshold"]

    def solidity(self, r_mean: float) -> float:
        pitch = 2.0 * math.pi * r_mean / self.n_vanes
        return self.chord / pitch

    def exit_flow_angle_deg(self, r_mean: float) -> float:
        """TE metal angle + Carter deviation. The flow does NOT leave at the metal angle."""
        camber = self.beta_le_deg - self.beta_te_deg
        return self.beta_te_deg + carters_deviation(camber, self.solidity(r_mean))

    def diffusion_factor(
        self, inlet: FlowState, V_out: float, Vt_out: float, r_mean: float
    ) -> float:
        """Lieblein's diffusion factor for a cascade.

            Df = 1 - V_out/V_in + (Vt_in - Vt_out) / (2 * sigma * V_in)

        Source: Lieblein, S., Schwenk, F.C. & Broderick, R.L. (1953), "Diffusion Factor
        for Estimating Losses and Limiting Blade Loadings in Axial-Flow-Compressor Blade
        Elements", NACA RM E53D01, Eq. 5. VERIFIED.

        VALIDITY: Df < 0.6. Above that the cascade SEPARATES, and every published loss
        correlation built on Df is fitted below it. RETRACTED (docs/centrifugal/
        17-tdd-plan.md slice S0): an earlier version of this docstring claimed that
        above Df >= 1 "the standard Lieblein wake-thickness correlation contains
        ln(1 - Df) and is UNDEFINED" -- that specific functional form could NOT be
        verified against the primary (NACA RM E53D01 is not held in
        research/papers/) and is now marked **UNVERIFIED**, not asserted. What IS
        established, independent of that functional form: Lieblein's own data, and
        every correlation descended from it, stop being fitted near Df ~ 0.6, and no
        published loss law exists above that range. See :meth:`loading_loss` for what
        this class does about running there anyway.
        """
        if V_out <= 0.0 or inlet.V <= 0.0:
            raise ValueError("diffusion factor needs positive velocities")
        sigma = self.solidity(r_mean)
        return 1.0 - V_out / inlet.V + (inlet.Vt - Vt_out) / (2.0 * sigma * inlet.V)

    def loading_loss(self, inlet: FlowState, V_out: float, Vt_out: float, r_mean: float) -> float:
        """Blade-loading (diffusion) loss: Δh_load = 0.05 * Df**2 * V_in**2.

        This is THE DOMINANT LOSS IN A HIGHLY-LOADED VANED DIFFUSER and its absence was the
        single largest error in the model: without it the whole diffusion system destroyed
        only 1.8% of P0 (5.4 kPa across the vanes) where it needs ~12%, and stage PR
        overshot NASA by 13%.

        Same correlation FORM as the impeller's :class:`ImpellerBladeLoadingCoppage`
        (Coppage 1956; Galvas, NASA TN D-7487 (1973), "Blade loading loss",
        Δh_BL = 0.05*D_f^2*u_3^2, verified verbatim), with the cascade's inlet velocity
        V_in in place of the rotor's blade speed U2 -- because a stationary row has no
        blade speed, and the loading is set by the velocity the flow actually arrives with.
        The coefficient 0.05 is UNCHANGED from Galvas.

        ⚠️ RANGE. At the HECC design point this diffuser runs at Df ~ 1.1, i.e. deep into
        separation and BEYOND THE RANGE OF EVERY PUBLISHED CASCADE CORRELATION. (An
        earlier version of this docstring additionally claimed "Lieblein's wake-
        thickness law is undefined at Df >= 1" -- that specific ln(1-Df) functional
        form is RETRACTED as UNVERIFIED, slice S0: NACA RM E53D01 is not held and the
        claim could not be checked against the primary. What stands regardless: every
        descendant correlation stops being fitted near Df ~ 0.6, well short of 1.1, so
        there is no published loss law here either way.) docs/PHYSICS-RULES.md: "A
        correlation outside its fitted range is an error, not an extrapolation. Guard
        and warn." So this warns -- loudly -- and reports the number rather than
        silently extrapolating. It does not pretend to be a validated prediction there.

        The loss is SELF-LIMITING, which is the physically right behaviour: more loss ->
        lower exit density -> higher exit velocity -> lower Df. The solve iterates to a
        consistent state rather than evaluating Df once on a loss-free guess.
        """
        Df = self.diffusion_factor(inlet, V_out, Vt_out, r_mean)
        if Df > self.df_warn_threshold:
            warnings.warn(
                f"cascade diffusion factor Df = {Df:.2f} exceeds Lieblein's separation "
                f"limit of 0.6 (V_in/V_out = {inlet.V / V_out:.2f}; real vaned diffusers "
                f"manage ~2.0-2.5). The blade-loading loss is being evaluated OUTSIDE the "
                f"range any published cascade correlation was fitted to"
                + (
                    " -- and above Df = 1 every descendant of Lieblein's correlation has "
                    "LONG SINCE stopped being fitted to anything (the fitted range ends "
                    "near Df ~ 0.6), so there is no published loss law at this condition "
                    "at all"
                    if Df >= 1.0
                    else ""
                )
                + ". The row is separated; treat this loss as an order-of-magnitude "
                "estimate, not a validated prediction.",
                RuntimeWarning,
                stacklevel=2,
            )
        return self.blade_loading_k * Df * Df * inlet.V * inlet.V

    def incidence_loss(self, inlet: FlowState) -> float:
        """Δh_inc = f_inc * (V * sin|alpha_in - beta_le|)**2 / 2.

        The velocity component normal to the vane's leading edge is taken as destroyed --
        the same model as the impeller inducer (Conrad), applied to a stationary row.
        """
        di = inlet.alpha - math.radians(self.beta_le_deg)
        W_star = inlet.V * math.sin(abs(di))
        return self.f_inc * 0.5 * W_star * W_star

    def friction_loss(self, inlet: FlowState, outlet_V: float, d_h: float) -> float:
        """Δh_sf = 2 * Cf * (L/d_h) * V_bar**2 -- the same form as Jansen's impeller
        skin friction, with the mean of inlet and exit velocity."""
        V_bar = 0.5 * (inlet.V + outlet_V)
        return 2.0 * self.Cf * (self.chord / d_h) * V_bar * V_bar


# --------------------------------------------------------------------- vaned diffuser


@dataclass(frozen=True)
class VanedDiffuser:
    """The vaned (channel) diffuser, r3 -> r4. A CASCADE -- never a rotor.

    Deviation by Carter's rule (:func:`carters_deviation`). ⚠️ A SLIP MODEL MUST NEVER BE
    APPLIED HERE (docs/PHYSICS-RULES.md rule 7): with omega = 0 the Coriolis term vanishes and the slip
    derivation collapses. This class accepts no slip model and there is no path to one.

    Losses (INTERNAL -- they destroy P0 and add no work, because a stationary row does no
    work):
      - incidence at the vane LE, off the DERIVED vane metal angle
      - skin friction along the vane passage
      - blade loading (Lieblein's diffusion factor) -- THE DOMINANT TERM

    Vane metal angles come from NASA Appendix C Tables C.25-C.28 (main + splitter vanes),
    extracted by my_scripts/extract_hecc_blade_angles.py. They are GEOMETRY, not knobs.

    ⚠️ STATION-DEPENDENT VANE COUNT (docs/centrifugal/17-tdd-plan.md slice S5;
    docs/centrifugal/20-review-fixes.md R10). HECC's diffuser splitter-vane leading edge
    sits at r = 0.2488 m, downstream of r3 = 0.23133 m -- only 20 (not 40) vanes exist
    over the first third of the passage. This class uses TWO different counts, each
    legal at its own station:
      - :attr:`_cascade` (n_vanes + n_splitters, 40) -- deviation ONLY. A
        TRAILING-EDGE quantity, and the splitters DO reach the TE.
      - :attr:`_cascade_loading` (:attr:`Z_eff_loading`, ~33.4575) -- the diffusion
        factor / loading loss (an LE -> TE quantity where the splitter is absent over
        the first third) AND skin-friction hydraulic diameter `d_h` (R10: `d_h` is a
        WETTED-PERIMETER quantity -- count the passages that actually exist over the
        chord, length-weighted exactly like the impeller's own `Z_eff` already does for
        the identical quantity class, data/coefficients.md S3.4a).
    See each property's docstring for why. Using 40 vanes everywhere (the pre-slice-S5
    behaviour) OVERSTATES solidity at the LE, which sits in Df's denominator, and so
    UNDERSTATES the loading loss. Using 40 (unweighted) for `d_h` (the pre-R10
    behaviour) likewise UNDERSTATES the friction loss relative to the impeller's own
    wetted-length convention.
    """

    r3: float
    r4: float
    b: float
    n_vanes: int
    beta_le_deg: float  # derived, Appendix C
    beta_te_deg: float  # derived, Appendix C
    n_splitters: int = 0
    blockage: float = 0.0
    Cf: float = FROZEN["diffuser.vaned.cf"]
    f_inc: float = FROZEN["diffuser.vaned.f_inc"]
    choke_k_low: float = FROZEN[
        "diffuser.choke.k_low"
    ]  # Aungier choke constants, reused at the diffuser throat
    choke_high_power: float = FROZEN["diffuser.choke.high_power"]
    choke_onset_scale: float = FROZEN["diffuser.choke.onset_scale"]
    choke_onset_threshold: float = FROZEN["diffuser.choke.onset_threshold"]
    # True vane chord, m, from NASA Appendix C (main vane 2.097 in = 0.053264 m). If None,
    # it is estimated from the radial extent and the mean metal angle -- a fallback, since
    # the real chord IS tabulated and should be used.
    chord: Optional[float] = None
    # SPLITTER VANE chord, m, from NASA Appendix C Tables C.27/C.28 (HECC: 1.411 in =
    # 0.0358396 m, my_scripts/extract_hecc_vane_angles.py). Feeds ONLY Z_eff_loading
    # below -- docs/centrifugal/17-tdd-plan.md slice S5. Ignored when n_splitters <= 0.
    chord_splitter: Optional[float] = None
    # SPLITTER VANE leading-edge radius, m (HECC: 0.248798 m, same source as chord_splitter).
    # Used ONLY as the fallback radial-extent proxy for Z_eff_loading when chord/
    # chord_splitter are not both supplied -- prefer the measured chords.
    splitter_le_r: Optional[float] = None
    # VANED-DIFFUSER THROAT AREA, m^2 -- the minimum passage between adjacent MAIN vanes.
    # THIS IS THE ELEMENT THAT CHOKES THE MACHINE. NASA's measured speedline chokes at
    # 5.24 kg/s; the inducer throat (0.020525 m^2) cannot be responsible -- it rotates, so
    # its capacity is set by relative-frame stagnation and is a hard, back-pressure-
    # independent ceiling at 5.92 kg/s. The diffuser's capacity is NOT fixed: it scales
    # with P02/sqrt(T02), and P02 collapses as the machine opens up. Evaluated at the P02
    # that exists AT choke it passes 5.18 kg/s against 5.24 measured.
    #
    # HECC: 0.005983 m^2 (9.27 in^2), 20 passages, from NASA Appendix C by minimum-distance
    # between adjacent vane surfaces (my_scripts/extract_hecc_throat_areas.py). The
    # SPLITTERS DO NOT BLOCK IT: the main->main throat cut makes ZERO crossings of the
    # splitter body (0.54 in clearance), and the two-piece cut through the splitter channels
    # is LONGER (0.928 vs 0.829 in) -- the HECC diffuser splitters begin downstream of the
    # throat. Assuming "20+20 = 40 passages" would HALVE the area and double the implied
    # choke margin, in the direction that looks like plausible blockage and gets tuned around.
    throat_area: Optional[float] = None
    name: str = "VanedDiffuser"

    @property
    def _resolved_chord(self) -> float:
        chord = self.chord
        if chord is None:
            chord = (self.r4 - self.r3) / math.cos(
                math.radians(0.5 * (self.beta_le_deg + self.beta_te_deg))
            )
        return chord

    @property
    def Z_eff_loading(self) -> float:
        """Length-weighted effective vane count for the DIFFUSION FACTOR / LOADING LOSS
        and for the SKIN-FRICTION HYDRAULIC DIAMETER `d_h` (docs/centrifugal/17-tdd-
        plan.md slice S5; docs/centrifugal/20-review-fixes.md R10) -- NOT for deviation
        (a TRAILING-EDGE quantity: the splitters DO reach the TE, so :attr:`_cascade`
        keeps the full ``n_vanes + n_splitters`` there, and ONLY there).

        WHY `d_h` USES THIS COUNT TOO (R10). A hydraulic diameter is a WETTED-PERIMETER
        quantity: it must count the passages that actually EXIST over the chord being
        integrated, length-weighted for a splitter that only occupies part of it --
        exactly this construct. Before R10, the friction `d_h` used the unweighted
        ``n_vanes + n_splitters = 40``, the OPPOSITE convention from the impeller's own
        skin-friction `d_h`, which already uses the fractional, length-weighted
        ``Impeller.Z_eff`` for the identical quantity class one component upstream
        (data/coefficients.md S3, row 189). Two conventions for the same physical
        quantity, in two components that sit next to each other, and both prior choices
        (diffuser: 40; impeller: already-correct Z_eff) happened to be the
        loss-INCREASING option -- see data/coefficients.md S3.4a for the measured effect
        of applying this one consistently.

        THE BUG THIS FIXES (the loading loss, slice S5). HECC's diffuser splitter-vane LEADING EDGE sits at
        r = 0.2488 m, while the vaned diffuser spans r3 = 0.23133 m -> r4 = 0.28459 m --
        over the FIRST THIRD of the passage only the 20 main vanes exist, exactly where
        the circulation per vane (and so the loading) is largest. Solidity sits in
        Lieblein's Df DENOMINATOR, so crediting 40 vanes there UNDERSTATES the loading
        loss.

        Same Aungier-derived formula as :attr:`Impeller.Z_eff`, one component
        downstream (Yang, Liu & Zhao 2023, *Machines* 11(1):118, Eq. (1) -- NOT "Li et
        al.", docs/centrifugal/17-tdd-plan.md slice S0):

            Z_eff,vd = n_vanes + n_splitters * (L_splitter / L_main)

        with L_main/L_splitter the vane CHORDS (LE -> TE), not a radial-extent proxy.
        Preferred: ``chord``/``chord_splitter`` from NASA Appendix C Tables C.25-C.28
        (my_scripts/extract_hecc_vane_angles.py). HECC: L_main = 0.0532632 m
        (2.097 in), L_splitter = 0.0358396 m (1.411 in) -> Z_eff,vd = 33.4575.

        Bit-identical to ``n_vanes`` when ``n_splitters <= 0`` (no splitter to weight) --
        Eckardt-scale unsplittered machines are untouched.

        Fallback (chords not both supplied): the radial extent LE(splitter) -> TE over
        LE(main) -> TE, using ``splitter_le_r``. This is a PROXY, not the measured
        quantity -- callers who care must supply the chords (same discipline as
        ``Impeller.Z_eff``'s own docstring).
        """
        if self.n_splitters <= 0:
            return float(self.n_vanes)
        if self.chord is not None and self.chord_splitter is not None:
            if self.chord <= 0.0:
                raise ValueError("chord must be > 0")
            return self.n_vanes + self.n_splitters * (self.chord_splitter / self.chord)
        if self.splitter_le_r is not None:
            L_main = self.r4 - self.r3
            L_split = self.r4 - self.splitter_le_r
            if L_main > 0.0:
                return self.n_vanes + self.n_splitters * (L_split / L_main)
        # No splitter geometry supplied at all -- the OLD (wrong) behaviour, preserved
        # only as a last resort so an under-specified caller does not crash.
        return float(self.n_vanes + self.n_splitters)

    @property
    def _cascade(self) -> _Cascade:
        """DEVIATION view: n_vanes + n_splitters (40 for HECC). A TRAILING-EDGE
        quantity and the splitters DO reach the TE -- legal, and unchanged by slice S5
        or R10.

        NOT used for friction's `d_h` since R10 (docs/centrifugal/20-review-fixes.md):
        `d_h` is a wetted-perimeter quantity and moved to :attr:`_cascade_loading`'s
        length-weighted count -- see :attr:`Z_eff_loading`'s docstring."""
        return _Cascade(
            beta_le_deg=self.beta_le_deg,
            beta_te_deg=self.beta_te_deg,
            n_vanes=self.n_vanes + self.n_splitters,
            chord=self._resolved_chord,
            Cf=self.Cf,
            f_inc=self.f_inc,
        )

    @property
    def _cascade_loading(self) -> _Cascade:
        """DIFFUSION FACTOR / LOADING LOSS view, AND (since R10) the FRICTION `d_h`
        view: :attr:`Z_eff_loading` (33.4575 for HECC, NOT 40). For the loading loss
        this is an LE->TE quantity, where the splitter is absent over the first third of
        the passage (docs/centrifugal/17-tdd-plan.md slice S5). For `d_h` it is the
        wetted-perimeter convention, matching the impeller's own friction `d_h`
        (docs/centrifugal/20-review-fixes.md R10)."""
        return _Cascade(
            beta_le_deg=self.beta_le_deg,
            beta_te_deg=self.beta_te_deg,
            n_vanes=self.Z_eff_loading,
            chord=self._resolved_chord,
            Cf=self.Cf,
            f_inc=self.f_inc,
        )

    def _choke_loss(self, inlet: FlowState, mdot: float, fluid: Air) -> float:
        """Aungier's choke-onset loss, applied to the DIFFUSER throat.

            Δh_ch = V_in**2 * (0.05*x + x**7) / 2,   x = 10*(1.1 - A_th/A*_th),  x > 0

        Same correlation and same published coefficients (0.05, the 7th power, the onset
        scale 10, the threshold 1.1) as :class:`ImpellerChokeAungier` (Kovář et al. 2021
        Eqs. (19)-(21), verified transcription of Aungier 2000) -- applied here, where the
        machine ACTUALLY chokes, instead of only at the inducer.

        ``A*_th`` is the standard 1-D critical (sonic) area for this mass flow at the
        diffuser inlet stagnation state:

            A*_th = mdot / (rho0 * a0) * [(gamma+1)/2] ** ((gamma+1)/(2*(gamma-1)))

        The loss switches on when the available throat falls to ~1.1x the area the flow
        needs to pass sonically, and rises as the 7th power thereafter -- which is what
        turns a flat speedline into the observed collapse. Returns 0 (and does NOT raise)
        when the throat is comfortably open: this is an onset loss, not a barrier.
        """
        if self.throat_area is None:
            return 0.0  # no throat supplied -> inert, exactly as before (debt D5)
        g = fluid.gamma
        rho0 = inlet.P0 / (fluid.R * inlet.T0)
        a0 = fluid.speed_of_sound(inlet.T0)
        A_star = (mdot / (rho0 * a0)) * ((g + 1.0) / 2.0) ** ((g + 1.0) / (2.0 * (g - 1.0)))
        x = self.choke_onset_scale * (self.choke_onset_threshold - self.throat_area / A_star)
        if x <= 0.0:
            return 0.0
        # V SQUARED -- see ImpellerChokeAungier. The published transcription omits the
        # square, which is dimensionally impossible (m/s, not J/kg).
        return inlet.V**2 * (self.choke_k_low * x + x**self.choke_high_power) / 2.0

    def solve(self, inlet: FlowState, mdot: float, fluid: Air) -> FlowState:
        # TWO cascade views, deliberately (docs/centrifugal/17-tdd-plan.md slice S5;
        # docs/centrifugal/20-review-fixes.md R10):
        #   cas       -- n_vanes + n_splitters (40 for HECC). A TRAILING-EDGE quantity:
        #                the splitters DO reach the TE, so this is legal for deviation
        #                (exit_flow_angle_deg) and incidence (keyed off the LE metal
        #                angle, not vane count).
        #   cas_load  -- Z_eff_loading (33.4575 for HECC, NOT 40). Used for (1) the
        #                diffusion factor and its loading loss, which span LE -> TE
        #                where the splitter is ABSENT over the first third of the
        #                passage (splitter LE at r=0.2488 m vs r3=0.23133 m ->
        #                r4=0.28459 m) -- crediting 40 vanes there understates the
        #                loading loss (solidity sits in Df's denominator) -- and (2),
        #                since R10, the skin-friction hydraulic diameter d_h below: d_h
        #                is a WETTED-PERIMETER quantity, so it counts the passages that
        #                actually exist over the chord, length-weighted exactly like the
        #                impeller's own friction d_h already does with Z_eff (both are
        #                the same Aungier construct one component apart --
        #                data/coefficients.md S3.4a).
        cas = self._cascade
        cas_load = self._cascade_loading
        r_mean = 0.5 * (self.r3 + self.r4)
        T0 = inlet.T0  # stationary row: T0 conserved, no work

        alpha4 = math.radians(cas.exit_flow_angle_deg(r_mean))
        area4 = 2.0 * math.pi * self.r4 * self.b * (1.0 - self.blockage)

        dh_inc = cas.incidence_loss(inlet)

        # Iterate: exit velocity depends on P0 loss, which depends on exit velocity.
        # d_h uses cas_load.n_vanes -- see the cascade-view note above (R10).
        d_h = (
            2.0
            * self.b
            * (2.0 * math.pi * r_mean / cas_load.n_vanes)
            / (self.b + 2.0 * math.pi * r_mean / cas_load.n_vanes)
        )
        # Iterate to a CONSISTENT state: the loading loss depends on the exit velocity,
        # which depends on the exit density, which depends on the loss. Evaluating Df once
        # on a loss-free guess would badly overstate it (the loss is self-limiting: more
        # loss -> lower rho -> higher V_out -> lower Df).
        P0 = inlet.P0
        for _ in range(200):
            Vm4 = _vm_at_angle(P0, T0, alpha4, area4, mdot, fluid)
            V4 = Vm4 / math.cos(alpha4)
            Vt4 = Vm4 * math.tan(alpha4)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # warn once, after convergence
                dh_load = cas_load.loading_loss(inlet, V4, Vt4, r_mean)
            dh_total = (
                dh_inc
                + cas.friction_loss(inlet, V4, d_h)
                + dh_load
                + self._choke_loss(inlet, mdot, fluid)
            )
            T4 = T0 - V4 * V4 / (2.0 * fluid.cp)
            P0_new = inlet.P0 * math.exp(-(dh_total / T4) / fluid.R)
            P0_new = P0 + 0.5 * (P0_new - P0)  # under-relax: the Df feedback is stiff
            if abs(P0_new - P0) < 1e-8 * P0:
                P0 = P0_new
                break
            P0 = P0_new

        # Re-evaluate once, unsuppressed, so the range warning fires on the CONVERGED Df.
        Vm4 = _vm_at_angle(P0, T0, alpha4, area4, mdot, fluid)
        V4 = Vm4 / math.cos(alpha4)
        cas_load.loading_loss(inlet, V4, Vm4 * math.tan(alpha4), r_mean)

        Vm4 = _vm_at_angle(P0, T0, alpha4, area4, mdot, fluid)
        Vt4 = Vm4 * math.tan(alpha4)
        out = state_from_totals(P0, T0, Vm4, Vt4, self.r4, self.b, fluid, s=0.0)
        return _with_entropy(out, inlet, fluid)


def _vm_at_angle(P0: float, T0: float, alpha: float, area: float, mdot: float, fluid: Air) -> float:
    """Vm satisfying continuity when the flow angle is FIXED (a bladed row sets the angle).

    Unlike the vaneless case, Vt is not independent here: Vt = Vm*tan(alpha).
    """

    def flux(Vm: float) -> float:
        V2 = Vm * Vm * (1.0 + math.tan(alpha) ** 2)
        T = T0 - V2 / (2.0 * fluid.cp)
        if T <= 1.0:
            return -1.0
        P = P0 * (T / T0) ** (fluid.gamma / (fluid.gamma - 1.0))
        return (P / (fluid.R * T)) * Vm * area - mdot

    lo, hi = 1e-6, 1.0
    for _ in range(200):
        V2 = hi * hi * (1.0 + math.tan(alpha) ** 2)
        T = T0 - V2 / (2.0 * fluid.cp)
        if T <= 1.0 or (hi >= fluid.speed_of_sound(max(T, 1.0))):
            break
        hi *= 1.05
    if flux(lo) * flux(hi) > 0.0:
        raise ValueError(
            f"no subsonic solution: mdot={mdot:.4f} kg/s, area={area:.6f} m2, "
            f"alpha={math.degrees(alpha):.1f} deg -- the row is CHOKED at this state"
        )
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if flux(lo) * flux(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------- exit guide vanes


@dataclass(frozen=True)
class ExitGuideVanes:
    """The deswirl row (EGV). A CASCADE -- Carter's rule, never slip.

    Turns the flow back toward axial. Metal angles from NASA Appendix C Tables C.29-C.30:
    LE 48.3 deg, TE 3.2 deg from meridional.

    ⚠️ DO NOT VALIDATE THIS AGAINST 34.3 deg. That number (data/hecc/design_point.csv,
    alpha_exit_deg) is the MEASURED rig swirl. NASA's DESIGN INTENT is 15 deg
    (docs/centrifugal/03-validation-matrix.md: "measured 34.3 deg against a design intent
    of 14 deg"). The rig UNDER-TURNS by ~20 deg relative to intent -- that is a real
    physical shortfall of the hardware, not a target for a geometry-driven model to hit.

    A model fed the true metal angles and a normal cascade deviation SHOULD land near
    15 deg. If it were "corrected" until it produced 34.3 deg, it would be reproducing a
    rig defect by mis-tuning the deviation rule -- fitting a model to a fault.
    """

    r_mean: float
    b: float
    n_vanes: int
    chord: float
    beta_le_deg: float
    beta_te_deg: float
    blockage: float = 0.0
    Cf: float = FROZEN["diffuser.egv.cf"]
    f_inc: float = FROZEN["diffuser.egv.f_inc"]
    name: str = "ExitGuideVanes"

    def solve(self, inlet: FlowState, mdot: float, fluid: Air) -> FlowState:
        cas = _Cascade(
            beta_le_deg=self.beta_le_deg,
            beta_te_deg=self.beta_te_deg,
            n_vanes=self.n_vanes,
            chord=self.chord,
            Cf=self.Cf,
            f_inc=self.f_inc,
        )
        T0 = inlet.T0
        alpha_out = math.radians(cas.exit_flow_angle_deg(self.r_mean))
        area = 2.0 * math.pi * self.r_mean * self.b * (1.0 - self.blockage)

        dh_inc = cas.incidence_loss(inlet)
        pitch = 2.0 * math.pi * self.r_mean / self.n_vanes
        d_h = 2.0 * self.b * pitch / (self.b + pitch)

        P0 = inlet.P0
        for _ in range(200):
            Vm = _vm_at_angle(P0, T0, alpha_out, area, mdot, fluid)
            V = Vm / math.cos(alpha_out)
            Vt_out = Vm * math.tan(alpha_out)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                dh_load = cas.loading_loss(inlet, V, Vt_out, self.r_mean)
            dh_total = dh_inc + cas.friction_loss(inlet, V, d_h) + dh_load
            T = T0 - V * V / (2.0 * fluid.cp)
            P0_new = inlet.P0 * math.exp(-(dh_total / T) / fluid.R)
            P0_new = P0 + 0.5 * (P0_new - P0)
            if abs(P0_new - P0) < 1e-8 * P0:
                P0 = P0_new
                break
            P0 = P0_new

        Vm = _vm_at_angle(P0, T0, alpha_out, area, mdot, fluid)
        Vt = Vm * math.tan(alpha_out)
        out = state_from_totals(P0, T0, Vm, Vt, self.r_mean, self.b, fluid, s=0.0)
        return _with_entropy(out, inlet, fluid)


__all__ = [
    "FlowState",
    "StationaryComponent",
    "VanelessSpace",
    "VanedDiffuser",
    "ExitGuideVanes",
    "carters_deviation",
    "state_from_totals",
    "solve_vm_from_continuity",
]
