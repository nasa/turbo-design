"""THE FROZEN COEFFICIENT SET — the single place every empirical constant lives.

**One coefficient set. All machines. No per-machine retuning.** (``docs/PHYSICS-RULES.md``,
validation doctrine.) A meanline code with ~10 loss correlations plus a slip model has enough
free parameters to hit *any single design point* by tuning. The freeze — not the fit — is what
makes a multi-machine result mean anything.

WHAT IS IN HERE. Every number that comes from a **correlation**: a loss coefficient, a slip
constant, a deviation rule, an empirical threshold that gates a correlation, and the gas
constants. Each one has a row in ``data/coefficients.md`` giving its source and its status
(VERBATIM / DERIVED / UNVERIFIED / PROXY).

WHAT IS **NOT** IN HERE, and must never be added:

* **Per-machine geometry** — ``r2``, ``b2``, ``backsweep_deg``, ``blockage``, ``tip_clearance``,
  ``Z_eff``, ``t1``, ``dbeta_dm_te``, throat area, blade counts. These vary by machine **by
  construction**. They belong on :class:`Impeller` / :class:`VanedDiffuser`, populated from the
  machine's own published coordinate tables. Freezing them would be a category error; *tuning*
  them per machine is not retuning, it is reading the geometry.
* **Numerical tolerances** — ``brentq``'s ``xtol``, iteration counts, bracket bounds. They set
  how precisely the solver answers, not what the answer is. ``data/coefficients.md`` §7.
* **Structural numbers** — the 2 in ``2*pi*r``, the 0.5 in an arithmetic mean, the 4 in ``D/4``.
  These are algebra, not physics, and citing them would be noise.

HOW THE FREEZE IS ENFORCED. ``tests/validation/test_frozen_coefficients.py`` asserts, in **both**
directions, that this dict and the ledger in ``data/coefficients.md`` §10 agree byte-for-byte:

* a value here that disagrees with the ledger **fails**;
* a key here that the ledger does not document **fails** (no undocumented knobs);
* a key the ledger documents that does not exist here **fails** (no phantom ledger rows).

The test contains **no literal coefficient values of its own**. It parses the ledger and imports
this module, and compares the two. (A regression guard that *copies* the configuration it guards
cannot fail for the reason it exists: this project once found one passing at ``rel=1e-12`` while
guarding a machine that no longer existed. ``docs/paper/00-INDEX.md`` §6.3.)

**After the freeze commit, changing any value below is a CI failure by construction, and it
invalidates the pre-registration in ``docs/centrifugal/21-prereg.md``.** Failures are REPORTED,
never tuned away. A coefficient changed to rescue a case has *hidden* something, not fixed it.

PROVENANCE OF THIS FILE. It was written by auditing every numeric literal in the module, and the
audit is itself a finding: **eleven of the coefficients below had no row in the ledger at all** —
the whole vaned-diffuser / EGV cascade (its friction coefficient, its incidence coefficient, and
**Carter's ``m = 0.26``, which sets the diffuser exit flow angle**), the diffuser's own reuse of
Aungier's choke constants at a different station, and the real-gas ``cp(T)`` table and Sutherland
constants that produce the polytropic efficiency this project reports against NASA. Their values
are UNCHANGED here; only their visibility is. ``data/coefficients.md`` §3.5 and §6b.
"""

from __future__ import annotations

from types import MappingProxyType
from typing import Mapping, Tuple, Union

Coefficient = Union[float, Tuple[Tuple[float, float], ...]]

FROZEN: Mapping[str, Coefficient] = MappingProxyType(
    {
        # ---------------------------------------------------------------- impeller: internal
        # Skin friction. Cf: Galvas, NASA TN D-7487 (1973). VERBATIM.
        # (The FORM of the loss is Oh, Yoon & Chung (1997) Table 6 -- Jansen (1967) is NOT
        # HELD and is not our authority for it. See data/coefficients.md, skin-friction row.)
        "impeller.skin_friction.cf": 0.004,
        # Conrad (1979) incidence. Kovář et al. 2021 Eq. (13), "commonly 0.6". VERBATIM.
        "impeller.incidence.f_inc": 0.6,
        # Coppage/Galvas blade loading, dh_BL = 0.05*Df^2*U2^2. Galvas TN D-7487. VERBATIM.
        "impeller.blade_loading.k": 0.05,
        # Galvas TN D-7487 Eq. (B59): K_BL = 0.75 conventional, 0.6 with splitters. VERBATIM.
        "impeller.diffusion_factor.k_bl_unsplittered": 0.75,
        "impeller.diffusion_factor.k_bl_splittered": 0.6,
        # Jansen (1967) clearance; Kovář et al. 2021 Eq. (37). VERBATIM (four confirmations).
        "impeller.clearance.jansen_coeff": 0.6,
        # Aungier (2000) mixing via Kovář et al. 2021 Eqs. (43)-(45); Lieblein D_eq stall threshold.
        "impeller.mixing.k": 0.5,
        "impeller.mixing.lieblein_deq_threshold": 2.0,
        # Aungier (2000) choke via Kovář et al. 2021 Eqs. (19)-(21). Ratio inverted vs the catalog's
        # literal transcription -- see data/coefficients.md §3.1.
        "impeller.choke.k_low": 0.05,
        "impeller.choke.high_power": 7.0,
        "impeller.choke.onset_scale": 10.0,
        "impeller.choke.onset_threshold": 1.1,
        # Aungier (2000) entrance diffusion via Kovář et al. 2021 Eqs. (22)-(24). Equation numbers
        # attributed to Meroni 2018 are `# UNVERIFIED -- paper NOT HELD` (R8).
        "impeller.entrance_diffusion.k": 0.4,
        # Gambini & Vellini (2021) Ch. 6 axial-length closure. Equation number UNVERIFIED
        # (paywalled). FALLBACK ONLY -- never reached when the machine supplies l_blade_m.
        "impeller.axial_length.gambini_c0": 0.014,
        "impeller.axial_length.gambini_c1": 0.023,
        "impeller.axial_length.gambini_c2": 1.58,
        # ---------------------------------------------------------------- impeller: parasitic
        # Daily & Nece (1960) disc friction via Kovář et al. 2021 Eqs. (52)-(53).
        # The leading 0.25*(.../4) double factor is a possible transcription artefact in Kovář et al..
        # REPRODUCED EXACTLY, NOT "FIXED" -- per the no-tuning rule.
        "impeller.disc_friction.leading_factor": 0.25,
        "impeller.disc_friction.denominator": 4.0,
        "impeller.disc_friction.laminar_c": 3.7,
        "impeller.disc_friction.turbulent_c": 0.102,
        "impeller.disc_friction.g_exponent": 0.1,
        "impeller.disc_friction.re_exponent": 0.2,
        # UNVERIFIED / disputed: Kovář et al. say 3e4, radcomp/CIMdes say 3e5. HECC Re2 ~ 3e7, so
        # turbulent either way. Immaterial here; NOT resolved. Matters for a small machine.
        "impeller.disc_friction.re_transition": 30000.0,
        # UNVERIFIED FOR THIS MACHINE: no measured back-face cavity exists. Daily & Nece tested
        # G = 0.0127-0.217; 0.05 is the order-of-magnitude midpoint. Weak knob (exponent 0.1).
        "impeller.disc_friction.back_face_clearance_ratio": 0.05,
        # Oh, Yoon & Chung (1997) recirculation; Kovář et al. 2021 Eq. (49). VERBATIM.
        # 90% of ALL parasitic work at the HECC design point -- debt D7, documented, NOT fixed.
        "impeller.recirculation.k": 8e-05,
        "impeller.recirculation.sinh_arg_coeff": 3.5,
        "impeller.recirculation.alpha2_power": 3.0,
        # A TRIPWIRE, NOT PHYSICS. Oh 1997 states no validity range for alpha2. UNVERIFIED --
        # engineering judgement. Fires at the HECC design point (alpha2 = 76.3 deg).
        # Do not cite it as Oh's.
        "impeller.recirculation.alpha2_warn_deg": 70.0,
        # Aungier gap discharge coefficient, sqrt(2/3). ⚠ The closure it feeds
        # (ImpellerLeakageAungier.delta_h) is DIMENSIONALLY BROKEN -- returns m/s, not J/kg
        # (debt D6, strict xfail). The value is unchanged pending Aungier (2000)'s own primary.
        # DO NOT GUESS THE MISSING FACTOR.
        "impeller.leakage.discharge_coeff": 0.816,
        # ---------------------------------------------------------------- slip
        # Wiesner (1967) via Kovář et al. 2021 Eq. (9). The industry-standard fit to Busemann.
        # ⚠ The Z fed to it is KNOB #6 (Z_eff = 25.4048, an UNVERIFIED CHOICE) -- but Z is
        # per-machine GEOMETRY, not a coefficient, so it lives on Impeller, not here.
        "slip.wiesner.exponent": 0.7,
        # UNVERIFIED (debt D2): Wiesner (1967) is paywalled; only the SHAPE of the two-branch
        # behaviour is corroborated. Inert today (r1_over_r2 defaults to None).
        "slip.wiesner.limiting_radius_ratio_c": 8.16,
        "slip.wiesner.limiting_radius_ratio_power": 3.0,
        # Stanitz (1952). CITED-ONLY (consistent across every secondary source).
        "slip.stanitz.k": 0.63,
        # Qiu et al. (2011) J. Turbomach. 133:041018. The leading 2 is VERBATIM from the
        # APPENDIX's own derivation; the paper's body (Eq. 6) omits it -- a typo, not an
        # alternative reading. Eq. (10b)'s 4*cos(beta2b) denominator likewise.
        "slip.qiu.f_leading_factor": 2.0,
        "slip.qiu.turn_denominator": 4.0,
        # ---------------------------------------------------- diffusion system (diffusion.py)
        # ⚠ EVERY KEY IN THIS BLOCK WAS ABSENT FROM THE LEDGER UNTIL THE FREEZE AUDIT.
        # Values unchanged; only their visibility. data/coefficients.md §3.5.
        #
        # Carter, A.D.S. (1950), ARC CP 29; Dixon & Hall 7th ed. Eq. 3.31 / Fig. 3.16.
        # THE DEVIATION RULE FOR BOTH STATIONARY CASCADES -- it sets the diffuser and EGV EXIT
        # FLOW ANGLE directly, and it was the single most consequential undocumented constant
        # in the model. # UNVERIFIED AT THIS PRECISION: Carter's m is a function of stagger and
        # is read from a chart; 0.26 is the widely-quoted circular-arc value.
        # NOT a slip model (PHYSICS-RULES rule 7): at omega = 0 the Coriolis term vanishes.
        "diffuser.carter.m": 0.26,
        # Smooth-wall skin-friction coefficient. # UNVERIFIED against a primary (the docstring
        # cites Japikse 1996 Ch. 4, NOT HELD; same order as Aungier's 0.005 for a vaneless
        # space). NOTE it differs from the impeller's Jansen Cf = 0.004, and nothing sources
        # the difference. Four independently-settable fields carry it.
        "diffuser.vaneless.cf": 0.005,
        "diffuser.cascade.cf": 0.005,
        "diffuser.vaned.cf": 0.005,
        "diffuser.egv.cf": 0.005,
        # Incidence coefficient, "same basis as the impeller's Conrad f_inc" -- an inference,
        # not a citation. # UNVERIFIED for a stationary row.
        "diffuser.cascade.f_inc": 0.6,
        "diffuser.vaned.f_inc": 0.6,
        "diffuser.egv.f_inc": 0.6,
        # The SAME Galvas 0.05 as the impeller's blade loading, reused for a cascade with the
        # inlet velocity in place of U2. Coefficient unchanged from Galvas.
        "diffuser.cascade.blade_loading_k": 0.05,
        # Lieblein separation limit. A GUARD, not a term: above it the cascade is separated and
        # every published loss correlation built on Df was fitted below it. HECC's diffuser runs
        # at Df ~ 1.1 -- OUTSIDE the fitted range of every published cascade correlation, which
        # this warns about loudly rather than silently extrapolating (PHYSICS-RULES: "a
        # correlation outside its fitted range is an error, not an extrapolation").
        "diffuser.cascade.df_warn_threshold": 0.6,
        # Aungier's choke constants REUSED at the DIFFUSER throat -- a physically distinct
        # station from the inducer throat the impeller's copies act at. This reuse was
        # undocumented. It is load-bearing: the DIFFUSER is what actually chokes on HECC
        # (5.82 kg/s), while the inducer choke term is inert at every achievable flow (R7).
        "diffuser.choke.k_low": 0.05,
        "diffuser.choke.high_power": 7.0,
        "diffuser.choke.onset_scale": 10.0,
        "diffuser.choke.onset_threshold": 1.1,
        # ------------------------------------------------------------------- gas (state.py)
        # Standard air. gamma is DERIVED (cp/(cp-R)) and is never an input: supplying cp, R and
        # gamma over-determines the gas and manufactures entropy across a LOSSLESS impeller.
        "air.cp": 1005.0,
        "air.R": 287.0,
        # ISA @ 288.15 K (ISO 2533:1975; White, Viscous Fluid Flow, Table 1.3).
        "air.mu": 1.81e-05,
        # Sutherland (1893); air constants as tabulated in White, 3rd ed., Eq. 1-36 / Table 1-2.
        # VERIFIED. Valid 100-1900 K. Was absent from the ledger until the freeze audit.
        "air.sutherland.mu_ref": 1.716e-05,
        "air.sutherland.t_ref": 273.15,
        "air.sutherland.s": 110.4,
        "air.sutherland.exponent": 1.5,
        # cp(T) for dry air, J/(kg K). Cengel & Boles, Table A-2b. TABULATED DATA, NOT A FIT.
        # This table produces `stage_eta_poly_realgas` -- the ONLY efficiency this project
        # reports against NASA's measurement -- and it, too, was absent from the ledger.
        "air.cp_table": (
            (250.0, 1003.0),
            (300.0, 1005.0),
            (350.0, 1008.0),
            (400.0, 1013.0),
            (450.0, 1020.0),
            (500.0, 1029.0),
            (550.0, 1040.0),
            (600.0, 1051.0),
        ),
    }
)


CANDIDATES: Mapping[str, Coefficient] = MappingProxyType(
    {
        # -------------------------------------------------------- CANDIDATE SETS (NOT FROZEN)
        # These belong to `turbodesign/centrifugal/candidates.py` -- loss sets used ONLY in the
        # R2 loss-channel experiment (`docs/centrifugal/25-prereg-r2.md`). They are declared in
        # their OWN mapping, separate from FROZEN above, so that:
        #   (a) candidates.py never carries a bare, chain-of-custody-less numeric literal --
        #       `tests/validation/test_frozen_coefficients.py`'s AST literal sweep covers
        #       candidates.py exactly like the frozen physics modules; and
        #   (b) they are NEVER compared against data/coefficients.md's BEGIN/END FROZEN SET
        #       block or the ledger-parity guards (1-3), because they are explicitly NOT
        #       frozen: `25-prereg-r2.md` §4 requires "every candidate is a new, named set
        #       with its own ledger rows" -- separate from the control, never mistaken for it.
        # Ledger rows: `data/coefficients.md` §12, "CANDIDATE SETS" (NOT the frozen-set block).
        #
        # Coppage (1956) recirculation, `dh_rc = k*sqrt(tan(alpha2))*Df^2*U2^2`, as transcribed
        # IDENTICALLY by two independent HELD sources: Kovář et al. 2021 Eq. (46); Galvas NASA
        # TN D-7487 (1973) Eq. (B73). Coppage (1956) itself is NOT HELD.
        # E-R2a, `ImpellerRecirculationCoppage` / `OhCoppageLossSet`.
        "candidate.recirculation_coppage.k": 0.02,
        # UNVERIFIED -- THIS PROJECT'S OWN. The arithmetic mean of the two distinct blade-count
        # populations in the 3-machine validation set (Eckardt O Z_exit=20, no splitters; HECC
        # and CC3 Z_exit=30, both 15+15 splittered) -- chosen BEFORE looking at what it does to
        # the answer (`25-prereg-r2.md` §3: "Do NOT search for a coefficient that makes it fit").
        # E-R2b Step 2, `ImpellerMixingBounded` / `MixingBoundedLossSet` -- A DIAGNOSTIC PROBE,
        # NEVER A VALIDATED MODEL, A CORRELATION, OR A RECOMMENDATION.
        "candidate.mixing_bounded.z_reference": 25.0,
    }
)


def candidate(key: str) -> Coefficient:
    """Read a candidate-set coefficient -- see :data:`CANDIDATES`. NOT part of the frozen set.

    Candidate coefficients belong to `turbodesign/centrifugal/candidates.py`'s experimental loss
    sets (`docs/centrifugal/25-prereg-r2.md`) and are deliberately excluded from the ledger-parity
    guards that protect :data:`FROZEN` -- see the module docstring above :data:`CANDIDATES`.
    """
    try:
        return CANDIDATES[key]
    except KeyError:  # pragma: no cover - a typo here fails every test that touches candidates.py
        raise KeyError(
            f"{key!r} is not a declared candidate coefficient. Every empirical constant used by "
            f"turbodesign/centrifugal/candidates.py must be declared in CANDIDATES here AND "
            f"documented in data/coefficients.md §12 (the CANDIDATE SETS section, NOT the "
            f"BEGIN/END FROZEN SET block)."
        ) from None


def frozen(key: str) -> Coefficient:
    """Read a frozen coefficient, failing loudly on a typo in the key.

    ``FROZEN[key]`` would raise ``KeyError`` anyway; this exists so the error names the
    ledger, because a missing key means the code and ``data/coefficients.md`` have diverged.
    """
    try:
        return FROZEN[key]
    except KeyError:  # pragma: no cover - a typo here fails every test that touches the model
        raise KeyError(
            f"{key!r} is not a frozen coefficient. Every empirical constant in "
            f"turbodesign/centrifugal/ must be declared here AND documented in "
            f"data/coefficients.md §10. Adding one without a ledger row fails "
            f"tests/validation/test_frozen_coefficients.py by construction."
        ) from None
