"""Candidate loss sets for the R2 loss-channel experiment (``docs/centrifugal/25-prereg-r2.md``).

**THIS MODULE IS NOT THE FROZEN CONTROL.** ``turbodesign.centrifugal.losses.OhLossSet`` stays
byte-identical and is never edited or imported-and-modified here -- every class below is a NEW,
separately-named object with its own ledger row. See ``docs/centrifugal/25-prereg-r2.md`` S0 for
why a candidate set is a new object, not an edit to the frozen one.

Every empirical constant this module needs (Coppage's ``k``, the mixing probe's ``z_reference``)
is declared in ``turbodesign.centrifugal.coefficients.CANDIDATES`` -- a mapping kept deliberately
SEPARATE from ``FROZEN``, so it is covered by the freeze guard's AST literal sweep
(``tests/validation/test_frozen_coefficients.py``) without ever being compared against the frozen
ledger block. Chain-of-custody rows for both constants: ``data/coefficients.md`` §12.

Two candidates, two different epistemic statuses:

1. :class:`ImpellerRecirculationCoppage` / :class:`OhCoppageLossSet` -- PUBLISHED, ZERO invented
   coefficients. Coppage (1956)'s recirculation loss, as transcribed identically by two
   independent HELD sources (Kovář et al. 2021 Eq. (46); Galvas NASA TN D-7487 (1973) Eq. (B73)).
   Swaps in for ``ImpellerRecirculationOh`` -- the PARASITIC (work) channel -- and nothing else.

2. :class:`ImpellerMixingBounded` / :class:`MixingBoundedLossSet` -- **NOT LITERATURE. THIS
   PROJECT'S OWN INVENTED FORM.** A diagnostic probe against the INTERNAL (pressure) channel's
   ``1/Z`` mixing collapse (``docs/centrifugal/24-the-map.md`` S2a). ``# UNVERIFIED`` in the
   loudest terms this project has -- see the class docstring. If it fits a measurement, that is a
   FIT, not a validation (``25-prereg-r2.md`` S3).

Every class here REUSES the frozen module's own machinery (``_diffusion_factor``,
``_camberline_length``, ``ImpellerLossState``, the individual frozen loss classes) by IMPORT, not
by copy -- the arithmetic those functions already implement is not re-derived a second time here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import ClassVar, Literal, Sequence

from .coefficients import CANDIDATES, FROZEN
from .losses import (
    ImpellerBladeLoadingCoppage,
    ImpellerChokeAungier,
    ImpellerClearanceJansen,
    ImpellerDiscFrictionDaily,
    ImpellerEntranceDiffusionAungier,
    ImpellerIncidenceConrad,
    ImpellerLeakageAungier,
    ImpellerLossState,
    ImpellerMixingAungier,
    ImpellerRecirculationOh,
    ImpellerSkinFrictionJansen,
    LossModel,
    LossSet,
    _camberline_length,
    _diffusion_factor,
)

# --------------------------------------------------------------- E-R2a: the WORK channel


@dataclass(frozen=True)
class ImpellerRecirculationCoppage:
    """Coppage (1956) impeller recirculation loss.

    **PROVENANCE (R8 rule -- never assert an equation number for a paper we do not hold).**
    Coppage's own primary (Coppage, J.E., Dallenbach, F., Eichenberger, J.P., Hlavaka, G.E.,
    Knoernschild, E.M., and Vanke, N., *Study of Supersonic Radial Compressors for Refrigeration
    and Pressurization Systems*, AiResearch Mfg. Co., WADC report) is **NOT HELD**. This class is
    transcribed from **two independent HELD sources that agree exactly on the coefficient and the
    functional form**:

        Kovář, P.; Tater, A.; Mačák, P.; Vampola, T. (2021), *Energies* 14(24), 8545, Eq. (46), p. 9 of 22:
            "Coppage in [8] estimated recirculation loss by a functional dependence on the exit
            absolute flow angle alpha2 and the diffusion factor D_f."

        Galvas, M.R. (1973), NASA TN D-7487, Eq. (B73), p. 51 (Galvas's station "3" = impeller
        exit = Kovář et al.'s station "2", confirmed against Galvas's own station table, p. 39):
            explicitly called "a *modified form* of the equation of reference 5" -- reference 5
            being the same Coppage & Dallenbach report Kovář et al. cite -- so the two are independent
            transcriptions of the same underlying (unheld) source, not one copying the other.

    Both give the identical leading coefficient (0.02) and functional form:

        dh_rc = 0.02 * sqrt(tan(alpha2)) * Df^2 * U2^2,     alpha2 in RADIANS

    ``D_f`` is the SAME diffusion factor Oh's own recirculation term and the frozen
    :class:`~turbodesign.centrifugal.losses.ImpellerBladeLoadingCoppage` already use --
    :func:`turbodesign.centrifugal.losses._diffusion_factor`, IMPORTED here, not reimplemented
    (research notes 01 S4.2: "D_f ... drives #4 and #9").

    PARASITIC (docs/PHYSICS-RULES.md rule 6): backflow from the diffuser into the impeller tip
    that consumes shaft work with no static-pressure benefit -- same physical mechanism as
    ``ImpellerRecirculationOh``, a different (and, per ``research/notes/10-recirculation-models.md``,
    much milder-growing) functional form of the SAME thing.

    **Sub-exponential, against Oh's `sinh` of a cubic**
    (``research/notes/10-recirculation-models.md`` "Ranking" table): Coppage's `sqrt(tan alpha2)`
    diverges as `alpha2 -> 90 deg`, but as a power law, not an exponential -- and Kovář et al.'s own
    error-minimising search selected Coppage over Oh in ALL THREE of their best-fitting loss sets
    (``research/notes/10-recirculation-models.md`` S2).

    ⚠ **THE DOMAIN GUARD IS OURS, NOT SOURCED -- and it is LOUD, not silent.**
    ``sqrt(tan(alpha2))`` is undefined (complex) for ``alpha2 <= 0`` or ``alpha2 >= 90 deg``.
    **Neither Kovář et al. nor Galvas states a guard anywhere** (``research/notes/10-recirculation-
    models.md`` S2, "Neither Kovář et al. nor Galvas states a guard; this is an inference from the
    formula's own algebra, not a printed validity range"). This class RAISES ``ValueError`` --
    loudly, at the call site, with the offending angle printed -- rather than silently clamping.
    A silent clamp is exactly the "term that cannot be caught being wrong" defect this project has
    now shipped eight times (docs/PHYSICS-RULES.md). `# UNVERIFIED -- this project's own addition;
    neither held source specifies any guard at all, silent or otherwise.`
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "parasitic"
    k: float = CANDIDATES["candidate.recirculation_coppage.k"]  # Coppage (1956), as transcribed
    # identically by Kovář et al. 2021 Eq. (46) and Galvas 1973 Eq. (B73) -- see class docstring.
    # NOT in FROZEN: this is a candidate set's own coefficient, chain-of-custody in
    # turbodesign/centrifugal/coefficients.py's CANDIDATES mapping and data/coefficients.md §12
    # (25-prereg-r2.md S0: a candidate set has its OWN ledger rows, and this experiment does not
    # add one to the frozen file -- see docs/centrifugal/25-prereg-r2.md).

    def delta_h(self, state: ImpellerLossState) -> float:
        Df = _diffusion_factor(state)
        alpha2 = math.atan2(
            state.Vt2, state.Cm2
        )  # radians, same convention as ImpellerRecirculationOh
        alpha2_deg = math.degrees(alpha2)
        if not (0.0 < alpha2_deg < 90.0):
            raise ValueError(
                f"ImpellerRecirculationCoppage: alpha2 = {alpha2_deg:.2f} deg is outside "
                "(0, 90) deg, where sqrt(tan(alpha2)) is undefined (complex). This guard is "
                "THIS PROJECT'S OWN addition -- # UNVERIFIED, neither Kovář et al. 2021 Eq. (46) nor "
                "Galvas 1973 Eq. (B73) states any validity range or domain guard for this "
                "formula (research/notes/10-recirculation-models.md S2). Raising loudly rather "
                "than silently clamping, per docs/PHYSICS-RULES.md's 'guard and warn, never "
                "silent' rule."
            )
        dh = self.k * math.sqrt(math.tan(alpha2)) * Df**2 * state.U2**2
        return max(dh, 0.0)


@dataclass(frozen=True)
class OhCoppageLossSet(LossSet):
    """``OhLossSet`` with ONLY the recirculation model swapped: Oh -> Coppage.

    Every other model is the frozen class, IMPORTED (not copied): incidence (Conrad), blade
    loading (Coppage), skin friction (Jansen), clearance (Jansen), mixing (Aungier), choke
    (Aungier), entrance diffusion (Aungier), disc friction (Daily & Nece), leakage (Aungier).
    Only ``ImpellerRecirculationOh`` -> :class:`ImpellerRecirculationCoppage`.

    This is E-R2a, the WORK-channel experiment (``docs/centrifugal/25-prereg-r2.md`` S2):
    recirculation is PARASITIC, so by rule 6 this swap can move psi/eta but MUST NOT move PR by
    more than ~1% on any machine. If it does, the internal/parasitic separation in this codebase
    is broken -- a far more serious finding than which recirculation correlation is better.
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
            ImpellerRecirculationCoppage(),  # <-- the ONE swap, everything else frozen/imported
            ImpellerLeakageAungier(),
        )
    )


# --------------------------------------------------------- E-R2b step 2: the INVENTED probe


@dataclass(frozen=True)
class ImpellerMixingBounded:
    """
    ⚠️ ============================================================================ ⚠️
    ⚠️  UNVERIFIED -- THIS PROJECT'S OWN INVENTED FORM, NOT A PUBLISHED CORRELATION.  ⚠️
    ⚠️  A DIAGNOSTIC PROBE TO LOCALISE A DEFECT. IF IT FITS, THAT IS A FIT, NOT A     ⚠️
    ⚠️  VALIDATION. IT MAY NEVER BE PRESENTED AS A VALIDATED MODEL, A CORRELATION,   ⚠️
    ⚠️  OR A RECOMMENDATION (docs/centrifugal/25-prereg-r2.md S3).                   ⚠️
    ⚠️ ============================================================================ ⚠️

    **What this is a probe FOR.** The frozen ``ImpellerMixingAungier`` computes
    ``dW = 2*pi*D2*Cu2 / (Z*L_tilde)`` -- Aungier's blade-to-blade circulation term. More blades
    (Z) => smaller dW => less mixing loss. Across this project's 3-machine validation set, Z_exit
    goes from 20 (Eckardt O, no splitters) to 30 (HECC, CC3 -- both 15+15 splittered), and mixing
    loss collapses 6.6x (``docs/centrifugal/24-the-map.md`` S2a: "mixing falls 6.6x via Aungier's
    dW ... more blades => smaller dW => less mixing loss"). That collapse is indicted as the
    mechanism behind the PR error tracking blade count. **There is NO published single-zone
    alternative in the held literature** (Oh's literal choice -- Johnston & Dean 1966 -- needs a
    two-zone jet/wake split this code cannot supply; ``docs/centrifugal/10-two-zone.md`` verdict:
    "do not implement"). So this class exists ONLY to test whether removing the ``1/Z`` collapse
    changes the sign/size of the pressure miss -- it is not offered as a better model of mixing.

    **What was changed, stated plainly.** Aungier's form is kept EXACTLY -- every other term
    (``W_sep``, ``Deq``, the Lieblein D_eq=2 stall threshold, the blocked exit area A2, the
    ``0.5*(W_sep-W_out)^2`` coefficient) is bit-for-bit the frozen ``ImpellerMixingAungier``
    arithmetic, reused via :func:`turbodesign.centrifugal.losses._camberline_length`. The ONLY
    change: ``Z`` in ``dW`` is held at a single fixed reference value, :attr:`z_reference`, for
    EVERY machine -- so Z can no longer vary between machines and the blade-count-driven collapse
    this term is suspected of causing cannot occur, by construction.

    **The reference value, and why it was chosen BEFORE looking at what it does to the answer.**
    ``z_reference = 25.0`` -- the arithmetic MEAN of the two distinct blade-count populations that
    exist in this project's validation set: Eckardt O's Z_exit=20 (no splitters), and HECC's /
    CC3's shared Z_exit=30 (both 15+15 splittered). This is the single most "neutral" way to
    remove Z as a per-machine free variable without favouring either population: it is not the
    mean weighted by machine COUNT (which would give 26.67, tilted toward the 2-of-3 splittered
    machines), it is the mean of the two DISTINCT VALUES that appear. It was NOT searched over,
    swept, or adjusted against any PR/psi/eta outcome -- ``docs/centrifugal/25-prereg-r2.md`` S3's
    explicit instruction: "Do NOT search for a coefficient that makes it fit. Use ONE value, for
    ALL machines, chosen for a STATED structural reason (not for what it does to the answer)."
    ``mixing_k`` and ``deq_threshold`` are UNCHANGED from the frozen set (Aungier's own published
    0.5 and Lieblein's 2.0), reused directly from ``FROZEN`` -- not re-tuned, not re-sourced.

    PARASITIC/INTERNAL: this is a MIXING loss -- INTERNAL (destroys P02, adds no work), identical
    classification to the frozen ``ImpellerMixingAungier`` it is a variant of.
    """

    kind: ClassVar[Literal["internal", "parasitic"]] = "internal"
    mixing_k: float = FROZEN["impeller.mixing.k"]  # UNCHANGED from the frozen set -- see docstring
    deq_threshold: float = FROZEN[
        "impeller.mixing.lieblein_deq_threshold"
    ]  # UNCHANGED from the frozen set
    z_reference: float = CANDIDATES[
        "candidate.mixing_bounded.z_reference"
    ]  # UNVERIFIED -- OURS. Mean of {20, 30}. See class docstring; chain-of-custody in
    # coefficients.py's CANDIDATES mapping and data/coefficients.md §12.

    def delta_h(self, state: ImpellerLossState) -> float:
        r2, b2 = state.r2, state.b2
        Z = self.z_reference  # <-- the ONE change from ImpellerMixingAungier.delta_h
        D2 = 2.0 * r2
        Cu2 = state.Vt2
        W2 = state.W2
        W1xi = math.hypot(state.Vm1, state.U1)
        L_tilde = _camberline_length(state)  # SAME closure ImpellerMixingAungier uses, imported

        dW = 2.0 * math.pi * D2 * Cu2 / (Z * L_tilde)
        W_max = (W1xi + W2 + dW) / 2.0
        Deq = W_max / W2
        W_sep = W2 if Deq <= self.deq_threshold else W2 * Deq / self.deq_threshold

        A2 = 2.0 * math.pi * r2 * b2 * (1.0 - state.blockage2)
        Wu1xi = state.U1
        W_out = math.hypot(state.Cm2 * A2 / (math.pi * D2 * b2), Wu1xi)

        dh = self.mixing_k * (W_sep - W_out) ** 2
        return max(dh, 0.0)


@dataclass(frozen=True)
class MixingBoundedLossSet(LossSet):
    """``OhLossSet`` with ONLY the mixing model swapped: Aungier's true-Z form ->
    :class:`ImpellerMixingBounded` (fixed-Z probe). Recirculation stays ``ImpellerRecirculationOh``
    (the FROZEN model, untouched) -- this set isolates the INTERNAL/pressure channel exactly as
    ``OhCoppageLossSet`` isolates the PARASITIC/work channel; the two defects are tested one at a
    time (``docs/centrifugal/24-the-map.md`` S2a: "correlated, not causally linked").

    ⚠️ **CARRIES THE SAME INVENTED-FORM WARNING AS** :class:`ImpellerMixingBounded` **-- READ IT.**
    Every table this loss set appears in must repeat: "UNVERIFIED -- THIS PROJECT'S OWN INVENTED
    FORM, NOT A PUBLISHED CORRELATION. A DIAGNOSTIC PROBE TO LOCALISE A DEFECT. IF IT FITS, THAT
    IS A FIT, NOT A VALIDATION."
    """

    models: Sequence[LossModel] = field(
        default_factory=lambda: (
            ImpellerIncidenceConrad(),
            ImpellerBladeLoadingCoppage(),
            ImpellerSkinFrictionJansen(),
            ImpellerClearanceJansen(),
            ImpellerMixingBounded(),  # <-- the ONE swap, everything else frozen/imported
            ImpellerChokeAungier(),
            ImpellerEntranceDiffusionAungier(),
            ImpellerDiscFrictionDaily(),
            ImpellerRecirculationOh(),  # FROZEN, unchanged -- this set isolates the internal channel
            ImpellerLeakageAungier(),
        )
    )


__all__ = [
    "ImpellerMixingBounded",
    "ImpellerRecirculationCoppage",
    "MixingBoundedLossSet",
    "OhCoppageLossSet",
]
