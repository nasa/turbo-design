"""Thermodynamic + kinematic state -- what the flow is doing.

The rothalpy-consistent core (``docs/centrifugal/01-design.md`` S3.2, fixing defect A
of ``docs/centrifugal/00-root-cause-analysis.md``):

    I       = h + W**2/2 - U**2/2                       # conserved across the impeller
    T0R2    = (I + U2**2/2) / Cp
    P0R2_is = P0R1 * (T0R2/T0R1)**(gamma/(gamma-1))      # <-- never just P0R1

``ThermoState.P0R_is`` makes that ideal relative total pressure a first-class,
independently-inspectable field, specifically so a test can assert it is not equal to
``P0R1`` (the bug this whole slice exists to catch).

Velocity triangles are meridional-referenced throughout: ``Vx = Vm*cos(phi)``,
``Vr = Vm*sin(phi)``, so ``Vm**2 == Vx**2 + Vr**2`` by construction (never
``Vr = W*sin(phi)`` -- defect E).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import ClassVar

from .coefficients import FROZEN


@dataclass(frozen=True)
class Air:
    """Constant-property ideal gas.

    TWO independent constants, never three. ``cp`` and ``R`` are primary; ``gamma`` is
    DERIVED. Do not "restore" gamma as an input field.

    An ideal gas is fixed by any two of (cp, cv, gamma, R), bound by

        cp = gamma*R/(gamma-1)      <=>      gamma = cp/(cp - R)

    Supplying all three OVER-DETERMINES the gas. If they do not satisfy that identity,
    the model silently **manufactures entropy from nothing**: an analytically isentropic,
    loss-free impeller closure reports

        ds = ln(T2/T1) * (cp - gamma*R/(gamma-1))

    The Slice 1 spec asked for cp=1005, gamma=1.4, R=287 -- which are inconsistent
    (1.4*287/0.4 = 1004.5, not 1005) and produced ds = +0.185 J/(kg K) across a LOSSLESS
    impeller. ``test_entropy_does_not_decrease`` caught it; nothing else did.

    Deriving gamma from (cp, R) fixes it at the source and keeps the two identities that
    matter algebraically exact:

        cp/R == gamma/(gamma-1)     (the isentropic P-T exponent)
        a    == sqrt(gamma*R*T)     (speed of sound -- needed for the Slice 10 choke
                                     criterion, which is on the MERIDIONAL Mach number)

    True gamma for cp=1005, R=287 is 1.39972, not 1.4. The difference is negligible in
    magnitude (~0.02% on the speed of sound) but the *inconsistency* is not: it is the
    difference between a model that enforces its thermodynamic identities and one that
    merely hopes they hold.

    The seam stays narrow so a Cantera-backed, T-dependent-cp gas can replace this later
    -- any such replacement MUST preserve the identities above.
    """

    cp: float = FROZEN["air.cp"]
    R: float = FROZEN["air.R"]
    # Dynamic viscosity, Pa*s. ISA standard-atmosphere value for air at 288.15 K
    # (ISO 2533:1975 / Sutherland's law reference point -- e.g. White, "Viscous Fluid
    # Flow", Table 1.3). Held CONSTANT, consistent with this class's existing
    # constant-cp/constant-R simplification. Needed from slice 3/4 onward for the
    # Reynolds numbers in the Jansen skin-friction and Daily & Nece disc-friction
    # correlations (turbodesign.centrifugal.losses); slices 1/2 never reference it.
    mu: float = FROZEN["air.mu"]

    def viscosity(self, T: float) -> float:
        """Dynamic viscosity at temperature T, by Sutherland's law, Pa*s.

            mu(T) = mu_ref * (T/T_ref)**1.5 * (T_ref + S)/(T + S),    S = 110.4 K

        Source: Sutherland, W. (1893), "The Viscosity of Gases and Molecular Force",
        Phil. Mag. 5(36), 507-531. Air constants (mu_ref = 1.716e-5 Pa*s at T_ref =
        273.15 K, S = 110.4 K) as tabulated in White, "Viscous Fluid Flow", 3rd ed.,
        Eq. 1-36 / Table 1-2. VERIFIED. Valid 100-1900 K.

        WHY THIS EXISTS. The constant ``mu`` above is the 288.15 K value, and it was
        being used at EVERY station. But the impeller exit runs at T2 ~ 400 K, where air
        is ~26% more viscous. Holding mu at its inlet value therefore OVER-PREDICTS the
        exit Reynolds number by ~1.26x -- the SAME ORDER as the entire suppressed-inlet
        effect (1/delta = 1.34x) that this model was separately criticised for ignoring.

        The error hid inside a "consistent with the constant-cp simplification"
        justification that was never true: cp varies a few percent over this range;
        mu varies by 26%. Dimensional and range discipline applies to fluid properties
        too, not just to the loss correlations that consume them.

        Anchor: mu(288.15) = 1.789e-5 against the ISA value 1.81e-5 (1.2%; the residual
        is Sutherland's own fit error at the reference point, not a coding error).
        """
        S, T_ref, mu_ref = (
            FROZEN["air.sutherland.s"],
            FROZEN["air.sutherland.t_ref"],
            FROZEN["air.sutherland.mu_ref"],
        )
        return mu_ref * (T / T_ref) ** FROZEN["air.sutherland.exponent"] * (T_ref + S) / (T + S)

    # Air cp(T), J/(kg K). VERIFIED tabulated values -- Cengel & Boles, "Thermodynamics:
    # An Engineering Approach", Table A-2b (dry air, ideal gas). Not a fit.
    _CP_TABLE: ClassVar[tuple] = FROZEN["air.cp_table"]

    def cp_at(self, T: float) -> float:
        """cp at temperature T, J/(kg K), by linear interpolation of the table above.

        THE CONSTANT ``cp`` IS THE 288 K VALUE, AND THE MACHINE DOES NOT STAY AT 288 K.
        HECC compresses to T02 ~ 484 K, where air's cp is 1026 -- 2.1% higher. That is not
        a rounding error, because polytropic efficiency is DEFINED through the exponent:

            eta_poly = [(gamma-1)/gamma] * ln(PR) / ln(T02/T01)

        so an error in gamma lands in eta ONE-FOR-ONE, with no physics in between. The
        mass-averaged cp over HECC's actual compression is 1012.8 (gamma = 1.3954, not
        1.3997), which INFLATES the reported eta by 0.8 points -- roughly HALF of the
        efficiency deficit this project spent a day hunting for in the loss models.

        The constant-cp simplification was documented, and being documented made it feel
        safe. It was not: a documented approximation is still an ERROR when the measurement
        you compare against did not make it.
        """
        tab = self._CP_TABLE
        if T <= tab[0][0]:
            return tab[0][1]
        for (T0, c0), (T1, c1) in zip(tab, tab[1:]):
            if T0 <= T <= T1:
                return c0 + (c1 - c0) * (T - T0) / (T1 - T0)
        return tab[-1][1]

    def cp_mean(self, T_a: float, T_b: float, n: int = 64) -> float:
        """Mass-averaged cp over a compression from T_a to T_b -- the cp that belongs in
        the polytropic exponent. Reduces to cp_at(T) as T_b -> T_a."""
        if abs(T_b - T_a) < 1e-9:
            return self.cp_at(T_a)
        return sum(self.cp_at(T_a + (T_b - T_a) * i / n) for i in range(n + 1)) / (n + 1)

    def gamma_mean(self, T_a: float, T_b: float) -> float:
        """gamma from the mass-averaged cp. gamma = cp/(cp - R) -- the identity is preserved."""
        cp = self.cp_mean(T_a, T_b)
        return cp / (cp - self.R)

    @property
    def gamma(self) -> float:
        """gamma = cp/cv = cp/(cp - R). DERIVED -- never an independent input."""
        return self.cp / (self.cp - self.R)

    @property
    def cv(self) -> float:
        return self.cp - self.R

    def speed_of_sound(self, T: float) -> float:
        """a = sqrt(gamma*R*T)."""
        return math.sqrt(self.gamma * self.R * T)


@dataclass(frozen=True)
class ThermoState:
    """The thermodynamic + kinematic state of the flow at one station.

    Frame discipline (``docs/PHYSICS-RULES.md`` rule 2): ``P0``/``T0`` are absolute-frame totals;
    ``P0R``/``T0R`` are relative-frame totals; ``P0R_is`` is the *ideal* (rothalpy-
    consistent, possibly loss-free) relative total pressure -- distinct from ``P0R``
    the moment a loss model subtracts something from it (slice 2+).
    """

    P0: float  # absolute total pressure, Pa
    T0: float  # absolute total temperature, K
    P: float  # static pressure, Pa
    T: float  # static temperature, K
    P0R: float  # actual relative total pressure, Pa
    P0R_is: float  # ideal (rothalpy-consistent) relative total pressure, Pa
    T0R: float  # relative total temperature, K
    rothalpy: float  # I = h + W**2/2 - U**2/2, J/kg
    U: float  # blade speed, m/s
    Vm: float  # meridional velocity, m/s
    Vx: float  # axial component of Vm, m/s
    Vr: float  # radial component of Vm, m/s
    Vt: float  # absolute tangential velocity, m/s
    W: float  # relative velocity magnitude, m/s
    s: float  # specific entropy relative to the stage inlet, J/(kg K)


def entropy(T: float, P: float, T_ref: float, P_ref: float, fluid: Air) -> float:
    """s - s_ref for a constant-Cp ideal gas: s = Cp*ln(T/T_ref) - R*ln(P/P_ref).

    Frame-independent (a function of the local static state only), so it is valid to
    evaluate on either the absolute or relative static state -- they are the same
    static point. Standard ideal-gas entropy relation (e.g. Cengel & Boles,
    *Thermodynamics*, Eq. 7-34).
    """
    return fluid.cp * math.log(T / T_ref) - fluid.R * math.log(P / P_ref)


def rothalpy(T: float, W: float, U: float, fluid: Air) -> float:
    """I = h + W**2/2 - U**2/2, with h = Cp*T for a constant-Cp ideal gas.

    Conserved across a rotor (docs/PHYSICS-RULES.md rule 3) -- the ``-U**2/2`` term is identically
    zero only when U is constant, i.e. only for an axial machine.
    """
    return fluid.cp * T + 0.5 * W**2 - 0.5 * U**2
