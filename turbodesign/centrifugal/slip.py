"""Slip factor models -- the centrifugal analogue of axial deviation.

``docs/PHYSICS-RULES.md`` rule 7: slip is a VELOCITY DEFICIT, not an additive deviation angle. It
enters the impeller-exit velocity triangle as

    Vt2 = sigma*U2 + Cm2*tan(beta2b)          (beta2b < 0 for backsweep)

never as ``beta2,flow = beta2b + delta`` (upstream's ``DeviationBaseClass`` returns
degrees, which would force a lossy angle<->velocity round trip -- see
``turbodesign/deviation/``, which is axial-only and NOT reused here).

Angle convention: every ``beta2b_deg`` in this module is the BLADE angle at the
impeller exit measured from the MERIDIONAL direction, negative for backsweep (the
HECC/CC3 convention -- see ``components.Impeller``). The classical slip papers
(Stodola 1927, Busemann 1928, Wiesner 1967) measure from the TANGENTIAL instead
(``sin(beta_tangential) == cos(beta_meridional)``); every formula below is already
written in the meridional convention, so no caller needs to convert.

Every coefficient below carries its source (author, year, equation). Anything not
traceable to a primary source is marked ``# UNVERIFIED`` -- see ``docs/PHYSICS-RULES.md``
("no magic numbers") and ``research/notes/01-centrifugal-meanline-methodology.md`` S3,
which is where every equation here was read from and verified (or not) against a PDF
in ``research/papers/``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Protocol, runtime_checkable

from .coefficients import FROZEN


class PhysicsError(Exception):
    """Raised when a physically meaningless computation is requested.

    The one case this module guards: applying a rotor slip model (driven by the
    Coriolis-induced relative eddy) to a stationary row. ``docs/PHYSICS-RULES.md`` rule 7 /
    Qiu et al. 2011 S2 (quoted in ``research/notes/
    01-centrifugal-meanline-methodology.md`` S6): *"for a radial vaned diffuser, it is
    not correct to use any existing slip models developed for radial impellers ...
    which essentially models only the first [Coriolis] term."* Use Carter's rule
    there instead (``research/notes/02-loss-models-catalog.md`` S13c).
    """


@runtime_checkable
class SlipModel(Protocol):
    """The slip-model plugin seam (``docs/centrifugal/01-design.md`` S3.3).

    ``sigma`` returns the Whitfield slip factor, ``sigma = 1 - C_slip/U2`` (research
    notes S3.1 -- the convention that stays correct in the presence of backsweep;
    ``Cu2 = sigma*U2 - Cm2*tan(beta_tilde_2)`` there uses the OPPOSITE sign convention
    for backsweep to this module's ``beta2b_deg`` -- this module and ``Stage.solve``
    use ``Vt2 = sigma*U2 + Cm2*tan(beta2b)`` with ``beta2b_deg`` already signed
    negative for backsweep, so the signs work out consistently within
    ``turbodesign.centrifugal``; do not mix the two sign conventions).
    """

    def sigma(
        self,
        beta2b_deg: float,
        Z: float,
        Cm2_over_U2: Optional[float] = None,
        omega: Optional[float] = None,
        **kw: object,
    ) -> float: ...


def _guard_rotating(omega: Optional[float]) -> None:
    """Refuse to compute slip for a non-rotating row.

    With ``omega = 0`` the Coriolis term that drives every slip model vanishes and the
    derivation collapses (``PhysicsError`` docstring above). ``omega=None`` means "not
    stated" and is deliberately NOT guarded, so geometry-only unit tests that never
    pass ``omega`` (e.g. ``test_wiesner_slip_factor``) keep working; only an EXPLICIT
    ``omega=0`` -- a vaned diffuser actually being asked to slip -- is refused.
    """
    if omega is not None and omega == 0.0:
        raise PhysicsError(
            "slip requires a rotating row (omega != 0): the Coriolis term vanishes at "
            "omega=0 and the slip derivation collapses (Qiu et al. 2011 S2). Use "
            "Carter's rule for a stationary (vaned diffuser) row instead."
        )


@dataclass(frozen=True)
class WiesnerSlip:
    """Wiesner (1967) -- the industry-standard fit to Busemann's exact solution.

    Basic correlation (research notes S3.5, verified via Kovář et al. 2021 Eq. 9):

        sigma_B = 1 - sqrt(cos(beta2b)) / Z**0.7

    Optionally applies Wiesner's limiting-radius-ratio correction if ``r1_over_r2`` is
    supplied (``r1_over_r2 = None`` -- the default -- skips it and returns ``sigma_B``
    unmodified, which is what every test in this slice exercises):

        eps_lim = exp(-8.16*cos(beta2b)/Z)
        sigma   = sigma_B * (1 - ((r1_over_r2 - eps_lim)/(1 - eps_lim))**3)   if r1_over_r2 > eps_lim
                = sigma_B                                                     otherwise

    ⚠ UNVERIFIED: the constants ``8.16`` and the exponent ``3``. Wiesner (1967) is
    paywalled (ASME) and these could not be obtained from a primary source -- only the
    *existence and shape* of the two-branch behaviour is verified (ETC2023-167;
    Qiu 2011 S3.1). Do not trust the magnitude of the correction; verify against
    Wiesner (1967), Dixon & Hall S7, or Aungier (2000) before using it in anger.
    """

    exponent: float = FROZEN["slip.wiesner.exponent"]
    limiting_radius_ratio_c: float = FROZEN[
        "slip.wiesner.limiting_radius_ratio_c"
    ]  # UNVERIFIED: 8.16
    limiting_radius_ratio_power: float = FROZEN[
        "slip.wiesner.limiting_radius_ratio_power"
    ]  # UNVERIFIED

    def sigma(
        self,
        beta2b_deg: float,
        Z: float,
        Cm2_over_U2: Optional[float] = None,
        omega: Optional[float] = None,
        r1_over_r2: Optional[float] = None,
        **kw: object,
    ) -> float:
        _guard_rotating(omega)
        beta2b = math.radians(abs(beta2b_deg))
        sigma_b = 1.0 - math.sqrt(math.cos(beta2b)) / Z**self.exponent
        if r1_over_r2 is None:
            return sigma_b
        eps_lim = math.exp(-self.limiting_radius_ratio_c * math.cos(beta2b) / Z)  # UNVERIFIED: 8.16
        if r1_over_r2 <= eps_lim:
            return sigma_b
        return sigma_b * (
            1.0 - ((r1_over_r2 - eps_lim) / (1.0 - eps_lim)) ** self.limiting_radius_ratio_power
        )  # UNVERIFIED: exponent 3


@dataclass(frozen=True)
class StanitzSlip:
    """Stanitz (1952) -- derived for radial-bladed impellers, independent of beta2b.

        sigma = 1 - 0.63*pi/Z

    Source: research notes S3.3 (consistent across every secondary source and
    textbook checked; Stanitz 1952 itself was not obtained -- the 0.63 constant is
    CITED-ONLY, not independently re-derived here).
    Validity: ``-45deg <= beta2b <= +45deg`` approximately, best at ``beta2b ~ 0``.
    """

    k: float = FROZEN["slip.stanitz.k"]

    def sigma(
        self,
        beta2b_deg: float,
        Z: float,
        Cm2_over_U2: Optional[float] = None,
        omega: Optional[float] = None,
        **kw: object,
    ) -> float:
        _guard_rotating(omega)
        return 1.0 - self.k * math.pi / Z


@dataclass(frozen=True)
class BusemannSlip:
    """Busemann (1928) -- the "exact" 2-D inviscid conformal-mapping solution.

    Published only as charts (``sigma = sigma(beta2', Z, r1/r2)``), so it is not
    directly implementable from a formula. Research notes S3.4 (verified, ETC2023-167
    Introduction): *"The results of Busemann (1928) are treated as an 'exact' solution
    of the problem"*, and Wiesner (1967) is *"the industry standard fit"* to it -- the
    notes explicitly say to "use Wiesner (which is a fit to Busemann)" as the
    implementable stand-in. This class does exactly that: it delegates to Wiesner's
    BASIC correlation only (no radius-ratio correction -- that piece is specific to
    Wiesner's own two-branch extension, not part of the Busemann fit being cited).
    """

    exponent: float = FROZEN[
        "slip.wiesner.exponent"
    ]  # delegated copy of Wiesner's basic correlation

    def sigma(
        self,
        beta2b_deg: float,
        Z: float,
        Cm2_over_U2: Optional[float] = None,
        omega: Optional[float] = None,
        **kw: object,
    ) -> float:
        _guard_rotating(omega)
        beta2b = math.radians(abs(beta2b_deg))
        return 1.0 - math.sqrt(math.cos(beta2b)) / Z**self.exponent


@dataclass(frozen=True)
class QiuSlip:
    """Qiu, Japikse, Zhao & Anderson (2011) -- the unified, flow-dependent model.

    The DEFAULT slip model for this library: the only one of the four that responds to
    the operating point, via the exit flow coefficient ``phi2 = Cm2/U2`` (research
    notes S3.6 -- needed for off-design/speedline prediction, not just a design
    point). Full model (verified, research notes S3.6):

        sigma = 1 - dsigma_radial - dsigma_turn - dsigma_passage

        dsigma_radial  = F*pi*cos(beta2b)*sin(gamma2) / Z
        dsigma_turn    = F*s2*phi2/(4*cos(beta2b)) * (dbeta/dm)_2
        dsigma_passage = -F*phi2*s2*sin(beta2b)/(4*rho2*b2) * d(rho*b)/dm

        F = 1 - 2*sin(pi/Z)*sin(pi/Z + beta2b)*cos(beta2b)*sin(gamma2) - t2/(s2*cos(beta2b))

    RESOLVED (debt D3, docs/centrifugal/17-tdd-plan.md slice S0): the body of the paper
    (Eq. 6) prints ``F`` WITHOUT the leading ``2`` on the first product term; the
    Appendix prints it WITH the ``2``. This is NOT a coin-flip between two equally
    valid readings -- the Appendix's own DERIVATION produces the ``2`` step by step:

        OC = R2*[1 - 2*sin(pi/Z2)*sin(pi/Z2 + beta2b)*cos(beta2b)*sin(gamma2)]

    walking through that geometric construction term by term reproduces the ``2``
    exactly, which means the body's Eq. (6) is a **typo**, not an alternative
    derivation. **This class uses the ``2``**, and that choice is confirmed by the
    primary source's OWN derivation -- not merely by matching CIMdes's implementation
    (research notes S3.6), though it also does that. It additionally reproduces the
    paper's own sanity check (``gamma2=90deg, F=1, dbeta/dm=0`` reduces exactly to
    Stodola, ``sigma = 1 - pi*cos(beta2b)/Z``: at ``F=1`` the two forms only differ by
    the ``2``, and the check requires the coefficient that makes the *whole bracket*
    vanish at ``beta2b=0, Z->inf``, which needs the leading ``2``).

    Status: **VERBATIM (Appendix)**. Cite: Qiu, Japikse, Zhao & Anderson (2011),
    "A Unified Slip Factor Model for Axial and Radial Impellers", J. Turbomachinery
    133(4):041018, Appendix.

    ⚠ REDUCED IMPLEMENTATION, still. ``sigma()``'s core signature (matching
    ``SlipModel``) is ``(beta2b_deg, Z, Cm2_over_U2, omega)``; two more terms are
    carried through the protocol's ``**kw`` seam (``dbeta_dm``, ``r2`` -- see
    ``dsigma_turn`` below), and ``d(rho*b)/dm`` (passage-area gradient) is still not
    carried at all, because ``Impeller`` does not carry per-station passage-width
    geometry. Documented simplifications:

    1. ``gamma2 = 90deg`` (purely radial discharge, ``sin(gamma2) = 1``) -- true for
       HECC/CC3 and every case this library currently validates against; overridable
       via the ``gamma2_deg`` kwarg for a future mixed-flow impeller.
    2. ``t2/s2 = 0`` (zero blade thickness in ``F``) unless overridden via
       ``t2_over_s2`` -- omits the (typically small) thickness-blockage correction to
       the shape factor.
    3. ``dsigma_passage`` is dropped entirely (needs ``d(rho*b)/dm``, not available).
       HALF-RESOLVED (debt D1, slice S0): this is NOT a defect on top of the model's
       own reduced scope -- it is **Qiu's own validated practice**. Qiu 2011 S3.3
       states outright: *"it will always be assumed to be zero in all of our
       validation studies."* Qiu himself never carries this term either. For HECC the
       magnitude is negligible regardless (``|dsigma_passage| <= 0.006`` at the design
       point), but the point is not the size of the number -- it is that dropping it
       reproduces the primary source's methodology, not a simplification of it.
    4. ``dsigma_turn`` -- **RESOLVED, docs/centrifugal/17-tdd-plan.md slice S6.** Qiu
       2011 Eq. (10b), verbatim::

           dsigma_turn = F*s2*phi2*(dbeta/dm)_2 / (4*cos(beta2b))
           s2          = 2*pi*r2/Z     (Qiu's Nomenclature: "pitch at the blade exit")

       ``dbeta_dm`` (the SIGNED ``(dbeta/dm)_2``, rad/m, in THIS module's
       ``beta2b_deg`` sign convention -- negative for backsweep, see the module
       docstring) and ``r2`` (the exit radius, needed to form ``s2``) are accepted as
       kwargs, absorbed by every other model's ``**kw`` seam and therefore invisible
       to them. Absent ``dbeta_dm`` (the default, ``None``), ``dsigma_turn = 0`` --
       this is Qiu's OWN treatment of an unavailable gradient term (S3.3, quoted in
       point 3 above, for the sibling term ``dsigma_passage``): a term this model
       cannot evaluate is dropped, not proxied.

       **``Z`` in Qiu's OWN derivation is the EXIT blade count** (30 for HECC, both
       main and splitter rows -- ``Delta_theta = 2*pi/Z2`` is the angle between
       ADJACENT blades at the impeller EXIT, and at HECC's TE adjacent blades
       alternate main/splitter). It is the SAME ``Z`` this method's ``d_radial``
       term already uses; there is only one ``Z`` in this class's OWN formula, and
       it is an exit quantity throughout -- unlike ``_diffusion_factor``'s
       ``Z_exit`` vs. the inlet-only ``Z_le`` split in ``losses.py`` (a DIFFERENT
       station, a DIFFERENT loss model, not this one).

       **BUT: the ``Z`` this class is actually CALLED with is not 30.** R6
       (``docs/centrifugal/20-review-fixes.md``): the only production caller
       (``Stage.solve``, ``turbodesign/centrifugal/solver.py``) passes
       ``Z=Impeller.Z_eff`` -- the SAME length-weighted, wetted-area/hydraulic-
       diameter blade count that ``ImpellerSkinFrictionJansen`` uses and that
       ``WiesnerSlip`` (the shipped default) is itself called with, **not** the
       exit count 30 this docstring's formula derivation is stated in terms of
       (data/coefficients.md SS5, "knob #6"). For HECC, ``Z_eff = 25.4048``.

       This is a DELIBERATE, DOCUMENTED choice (option (a) of R6), not an
       oversight left unresolved: threading the physically-correct exit count
       (``Z=30``) into ``Stage.solve``'s ``self.slip.sigma(...)`` call would
       change ``WiesnerSlip``'s own sigma too (it is the SAME call site, the
       SAME ``Z`` kwarg, for whichever model is configured) -- from 0.9066 to
       0.9168 at the HECC design point, verified numerically to move psi and
       stage PR. ``WiesnerSlip`` is the shipped default in every validation
       fixture, so that is a real regression, not a latent one. ``Z_eff`` is
       therefore RETAINED for the slip call, exactly as it is for
       ``ImpellerSkinFrictionJansen``'s hydraulic diameter -- it is knob #6,
       reported, not turned.

       **Consequence for THIS class, since it is not yet the default:**
       ``sigma_qiu(Z=Z_eff=25.4048, real dbeta_dm=-4.8255) = 0.9265`` at the HECC
       design point -- NOT ``0.9361``, which is ``sigma_qiu`` at the physically-
       defensible ``Z=30`` (diagnostic only; data/coefficients.md SS5a). The
       18%-larger pitch ``s2 = 2*pi*r2/Z`` this produces is LATENT today only
       because Wiesner is the default slip model; it fires the moment ``QiuSlip``
       is substituted in (as it already is for the speedline candidate, slice 9).
       Both numbers are recorded, and neither is silently "the" value -- always
       state which ``Z`` a quoted ``sigma_qiu`` was computed at.

       **THE DELETED PROXY, for the record.** Before this slice, the unavailable
       ``dbeta/dm`` was stood in for by ``F*phi2*tan(beta2b)/4`` -- treating "total
       turning experienced" as the (positive-MAGNITUDE) exit blade angle itself. That
       proxy is not merely unverified: for HECC's S-shaped blade (Fig. 18: ``|beta|``
       has a mid-chord minimum and *rises* toward the TE) the TRUE signed
       ``(dbeta/dm)_2`` is **negative** (measured ≈ -4.8 rad/m at the RMS streamline,
       fillet-clear 65-88% chord window -- ``Impeller.dbeta_dm_te``'s own docstring;
       data/coefficients.md S5), while the magnitude-driven proxy has the WRONG SIGN
       and *inverts* how sigma responds to the exit flow coefficient (Qiu S3.2:
       negative ``dbeta/dm`` should make sigma *increase* with phi2; the proxy made it
       fall). Its ``sigma = 0.8854`` (HECC, as coded) landing near the "wanted" ≈0.89
       was a COMPENSATING ERROR, not agreement -- exactly why it is deleted outright,
       with no fallback, rather than patched.

    Sanity check preserved exactly: with ``Cm2_over_U2 = None`` (``phi2 -> 0``,
    ``dsigma_turn -> 0`` regardless of ``dbeta_dm``) and the defaults above, this
    reduces to ``sigma = 1 - F*pi*cos(beta2b)/Z`` -- Stodola scaled by the
    (thickness-free, radial-discharge) shape factor ``F``, matching the paper's own
    ``gamma2=pi/2, F=1 -> Stodola`` sanity check once ``F`` is further set to 1 (the
    ``Z -> infinity`` limit).

    Validity limit (Qiu 2011 S3.1, verified): ``r1/r2 < F``. Not checked here --
    ``Impeller`` does not carry ``r1`` in this slice's minimal signature.

    Wiesner remains the DEFAULT slip model for this library despite QiuSlip's greater
    physical content (docs/centrifugal/17-tdd-plan.md S6; data/coefficients.md S5):
    with the real, measured ``dbeta_dm`` and ``Z`` AS ACTUALLY WIRED by
    ``Stage.solve`` (``Z_eff = 25.4048``, R6 above), Qiu gives a HECC design-point
    sigma of **0.9265** -- further from the shipped model's headline (Wiesner,
    0.9066 at the same ``Z_eff``) than a naive read might suggest, and still
    *worse* agreement, not better: the fourth time an honest input has moved this
    model away from the measurement. (At the physically-defensible exit blade
    count ``Z=30`` -- not what is wired -- sigma rises further, to 0.9361;
    diagnostic only, data/coefficients.md S5a.) No model is adopted or rejected
    for what it does to the answer; QiuSlip remains the candidate for a
    SPEEDLINE (slice 9), where it is the only model that responds to the
    operating point at all.
    """

    f_leading_factor: float = FROZEN[
        "slip.qiu.f_leading_factor"
    ]  # Qiu et al. (2011) Appendix, verbatim
    turn_denominator: float = FROZEN["slip.qiu.turn_denominator"]  # Qiu et al. (2011) Eq. (10b)

    def sigma(
        self,
        beta2b_deg: float,
        Z: float,
        Cm2_over_U2: Optional[float] = None,
        omega: Optional[float] = None,
        gamma2_deg: float = 90.0,
        t2_over_s2: float = 0.0,
        dbeta_dm: Optional[float] = None,
        r2: Optional[float] = None,
        **kw: object,
    ) -> float:
        _guard_rotating(omega)
        beta2b = math.radians(abs(beta2b_deg))
        gamma2 = math.radians(gamma2_deg)
        phi2 = 0.0 if Cm2_over_U2 is None else Cm2_over_U2

        F = (
            1.0
            - self.f_leading_factor
            * math.sin(math.pi / Z)
            * math.sin(math.pi / Z + beta2b)
            * math.cos(beta2b)
            * math.sin(gamma2)
            - t2_over_s2 / math.cos(beta2b)
        )

        d_radial = F * math.pi * math.cos(beta2b) * math.sin(gamma2) / Z

        if dbeta_dm is None:
            # Qiu's own treatment of an unavailable gradient term (S3.3) -- NO
            # fallback, NOT the deleted tan(beta2b) proxy. See class docstring point 4.
            d_turn = 0.0
        else:
            if r2 is None:
                raise ValueError(
                    "QiuSlip.sigma: dbeta_dm was supplied without r2 -- Eq. (10b)'s "
                    "pitch s2 = 2*pi*r2/Z cannot be formed without the exit radius. "
                    "Supply both or neither."
                )
            s2 = 2.0 * math.pi * r2 / Z  # Qiu's Nomenclature: pitch at the blade exit
            d_turn = (
                F * s2 * phi2 * dbeta_dm / (self.turn_denominator * math.cos(beta2b))
            )  # Eq. (10b)

        return 1.0 - d_radial - d_turn
