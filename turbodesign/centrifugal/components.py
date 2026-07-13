"""Stage components -- currently just the impeller (slice 1 has no diffuser).

Per ``docs/centrifugal/01-design.md`` S3.5, a centrifugal component owns an LE and a
TE station (unlike upstream's single-cutting-line ``BladeRow``). For the impeller,
the LE station is the inducer eye (found from geometry, see
``turbodesign.centrifugal.geometry.MeridionalPath.inducer_eye``) and the TE station is
the requested exit radius.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Impeller:
    """A centrifugal impeller: blade count, backsweep, and exit radius.

    ``backsweep_deg`` is negative for a backswept blade (the HECC/CC3 convention;
    ``docs/centrifugal/00-root-cause-analysis.md`` S3.C), entering the no-slip exit
    triangle as ``Vt2 = U2 + Cm2*tan(backsweep_deg)``.

    ``blockage`` is the fraction of the *geometric* exit area assumed blocked by
    boundary layers/wakes, applied as ``area = area_geom*(1-blockage)``
    (``docs/PHYSICS-RULES.md`` rule 4). Default is **0.0**: slice 1 is pure geometry, with no
    empirical corrections. Real impellers DO have exit blockage (~0.9-0.95 effective
    area) from boundary layers and wakes, but that is an *empirical loss-model*
    quantity -- it must arrive in a later slice WITH A PRIMARY-SOURCE CITATION
    (``docs/PHYSICS-RULES.md``: "no magic numbers"), not as a silent geometric fudge here.

    A previous version of this module defaulted ``blockage=0.05``, justified as
    "back-solved to match the exit area NASA's coordinates imply (A2 = 0.020018)".
    That justification was circular: 2*pi*0.21581*0.015539 = 0.021071 m^2, and
    0.020018 = 0.95 * 0.021071 -- i.e. the cited "target" area already contains the 5%
    being solved for. The geometric area implied by NASA's own coordinates
    (0.021071 m^2) is already within 0.45% of NASA's analytic 2*pi*r2*b2 = 0.020977
    m^2, leaving no room for a further 5% blockage. Worse, it was a free knob no test
    could see: blockage 0.00/0.05/0.08 give PR 6.817/6.734/6.679, and all three pass
    the tests' +-2% band -- a fudge that would have silently absorbed slip error once
    slip lands (slice 2).

    ``x_le`` is the axial location of the blade leading edge (impeller inlet / inducer
    eye station). If ``None`` (the default), it falls back to
    :meth:`turbodesign.centrifugal.geometry.MeridionalPath.inducer_eye`'s
    argmin-hub-radius heuristic -- a MODELLING ASSUMPTION, not a geometric fact, that
    is benign in slice 1 but becomes load-bearing from slice 4 (incidence loss, inducer
    relative Mach) onward. Prefer supplying it explicitly once it is known.
    """

    n_blades: int
    backsweep_deg: float
    r_te: float
    blockage: float = 0.0
    x_le: Optional[float] = None
    n_splitters: int = 0
    splitter_le_r: Optional[float] = None
    # Radial (shroud) blade tip gap, m -- the unshrouded clearance ``eps`` feeding the
    # Jansen clearance loss (research notes 02 S6a) and the Aungier leakage loss
    # (research notes 02 S10a) in turbodesign.centrifugal.losses. Default 0.0
    # (shrouded / no-clearance impeller, matching slices 1-4's assumption). Slice 5's
    # HECC point supplies 0.000305 m = 0.012 in, NASA CR-2014-218114/REV1 Table 1
    # (verified).
    tip_clearance: float = 0.0
    # TRUE meridional camberline lengths (LE -> TE), m -- the lengths Aungier/Li's
    # splitter formula actually calls for. When supplied they are used directly and the
    # radial-extent proxy below is not consulted. For HECC these come from NASA's own
    # blade coordinates (Appendix C Tables C.3-C.13 main, C.14-C.24 splitter), extracted
    # by my_scripts/extract_hecc_blade_angles.py: span-averaged L_FB = 0.20988 m
    # (8.2628 in), L_SB = 0.14559 m (5.7315 in).
    l_main_m: Optional[float] = None
    l_splitter_m: Optional[float] = None
    # TRUE along-blade (camberline) length, m -- docs/centrifugal/17-tdd-plan.md slice S1.
    # A DIFFERENT quantity from l_main_m (the MERIDIONAL length, LE->TE in the x-r plane,
    # S6.4/S7.2 of the closure audit): l_blade_m is the length ALONG the twisted blade
    # surface, always >= l_main_m, since a blade that turns the flow is longer than its
    # meridional projection. It is what Jansen's skin friction acts over and what
    # Aungier's mixing loss's blade-to-blade Delta-W integrates over -- both wetted-
    # passage-length quantities, not meridional-only ones
    # (turbodesign.centrifugal.losses._camberline_length).
    #
    # ADOPTED for HECC: the MEASURED 3-D camberline arc length, ds = sqrt(dm^2 +
    # (r*dtheta)^2), integrated along each of NASA Appendix C's 11 main-blade sections
    # (Tables C.3-C.13, data/hecc/blade_angles.csv) and span-averaged
    # (my_scripts/extract_hecc_blade_angles.py's blade_length_estimates) = 0.237875 m.
    # It needs NO blade angle at all -- it is the geometry, directly. This is the
    # value my_scripts/hecc_stage.py's IMPELLER dict actually ships
    # (l_blade_m=0.237875), and it is guarded by
    # tests/unit/test_camberline_length.py::
    # test_blade_length_is_the_measured_arc_not_the_LE_angle_projection.
    #
    # REJECTED: l_main_m / mean(cos(beta1b_deg)), span-averaged over the same 11
    # sections = 0.209875 / 0.717782 = 0.292394 m. This divides the WHOLE blade's
    # meridional length by the cosine of the INDUCER LEADING-EDGE blade angle (45.46
    # deg) -- A ONE-STATION QUANTITY (beta1b_deg is measured only at the LE) APPLIED
    # TO THE WHOLE BLADE, exactly the cross-station-proxy error this project keeps
    # finding elsewhere (Z_le, A_th, Lm -- docs/PHYSICS-RULES.md). HECC's blade is
    # S-shaped (NASA Fig 18: a mid-chord beta minimum near 14-20 deg, where
    # cos(beta) ~ 0.95, not the LE's ~0.70), so this proxy inflates the length by
    # ~23% relative to the measured arc, and the skin-friction/mixing losses with it.
    #
    # This projection was very nearly adopted anyway, for exactly the wrong reason:
    # it happens to agree with a rough L_m/<cos beta> figure this project's own prior
    # research sessions had already used to validate the S1 prediction against
    # (~0.2925 m). THAT IS NOT VALIDATION. Agreeing with your own prior estimate is
    # not an independent check -- it is confirmation bias wearing a docstring. The
    # measured arc (0.237875 m) is adopted because it is geometry, derived once,
    # directly, from NASA's own coordinates; the LE-angle projection (0.292394 m) is
    # rejected because it is a derivation whose real selling point was that it
    # matched a number this project already believed, which is the
    # derivation-fitted-to-prediction anti-pattern, not a property of the geometry.
    #
    # IDENTIFICATION, marked honestly: Jansen (1967) itself is not held as a primary
    # source, so "L_blade is the wetted-passage length Jansen's Cf acts over" is an
    # INFERENCE from his friction form's structure (Δh_sf = 2*Cf*(L/d_h)*W_bar^2), not
    # a citation to Jansen's own definition. # UNVERIFIED against Jansen (1967) primary.
    l_blade_m: Optional[float] = None
    # Main-blade LE radius, m -- used ONLY by the Z_eff radial-extent fallback, so that a
    # blade is measured from its leading edge rather than from the shaft centerline. If
    # the true meridional lengths above are supplied, this is not consulted. HECC hub LE:
    # 0.040485 m (1.5939 in), NASA Appendix C Table C.3.
    r_le_hint: Optional[float] = None
    # Inducer BLADE METAL ANGLE at the RMS inlet streamline, degrees from meridional,
    # positive in the same sense as the LE relative flow angle atan2(U1, Vm1).
    #
    # This is the field whose ABSENCE made the incidence loss structurally zero. With no
    # blade angle to compare the flow against, ImpellerIncidenceConrad had to define
    # beta_opt := beta_flow, which makes W* == 0 identically -- a loss that returns
    # EXACTLY 0.0 at every operating point, forever. That is the same silent-zero defect
    # (docs/PHYSICS-RULES.md defect D) this module exists to avoid, and it means efficiency cannot
    # roll over on a speedline: incidence IS the surge-side rollover.
    #
    # It is NOT tabulated by NASA (Table 2 gives only LEAN angles, a different quantity)
    # but it IS derivable from NASA's own blade coordinates -- see
    # my_scripts/extract_hecc_blade_angles.py and data/hecc/blade_angles.csv. HECC at the
    # RMS streamline (60.8% span): 45.46 deg.
    #
    # Evaluate it at the RMS streamline, because that is where the loss state's U1/Vm1
    # live. A spanwise distribution is not representable in a 1-D meanline model.
    inducer_blade_angle_deg: Optional[float] = None
    # LE blade TANGENTIAL thickness, m -- the ``t_1xi`` of Conrad's blockage-corrected
    # optimum-incidence angle. Default 0.0 collapses the blockage correction to the
    # identity (beta_opt == beta1b), which is the correct BOUNDARY CASE, not a fudge: in
    # NASA's tabulated sections the two blade surfaces meet at the LE, so the tabulated
    # LE thickness genuinely IS zero. A real blade has a finite (elliptical) nose; supply
    # it if known. Second-order next to the blade angle itself.
    le_blade_thickness: float = 0.0
    # INDUCER THROAT AREA, m^2 -- the minimum passage between adjacent MAIN blades near
    # the LE. Debt D5 closed: ImpellerChokeAungier used to stand A_th in as A1 (the LE
    # annulus), which OVERSTATES the real throat by 53% and made the choke loss
    # structurally inert -- the model could not choke at all, so its speedline stayed flat
    # where NASA's collapses 20 efficiency points.
    #
    # HECC: 0.020525 m^2 (31.81 in^2), from NASA Appendix C by minimum-distance between
    # adjacent blade surfaces (my_scripts/extract_hecc_throat_areas.py; cross-checked
    # against the textbook o = pitch*cos(beta1b) - t_LE to within 1-6% at every span).
    # 15 passages -- the SPLITTERS DO NOT BLOCK THIS THROAT (verified geometrically: the
    # throat cut clears the splitter by 1.4-2.6 in at every span). Assuming "15+15=30"
    # would halve the area and double-count blockage.
    throat_area: Optional[float] = None
    # SIGNED (dbeta/dm)_2, rad/m, at the impeller EXIT -- docs/centrifugal/
    # 17-tdd-plan.md slice S6, feeding QiuSlip's dsigma_turn (Qiu 2011 Eq. 10b). SAME
    # sign convention as backsweep_deg (negative for backsweep): Qiu S3.2 states that
    # a blade angle DECREASING toward the exit (signed dbeta/dm < 0, in this
    # convention) makes sigma INCREASE with the exit flow coefficient -- HECC's own
    # S-shaped blade (|beta| minimum mid-chord, RISING toward the TE, NASA Fig 18)
    # satisfies exactly that: signed beta becomes MORE negative toward the TE.
    #
    # Default None -- QiuSlip's own treatment of an absent gradient (Qiu S3.3's own
    # practice for the sibling dsigma_passage term): dsigma_turn = 0, NOT the deleted
    # tan(beta2b) proxy this slice removes (turbodesign.centrifugal.slip.QiuSlip).
    # WiesnerSlip (the DEFAULT slip model for every fixture) ignores this field
    # entirely via its **kw seam, so populating it is output-INVARIANT.
    #
    # HECC: -4.8255 rad/m (my_scripts/extract_hecc_blade_angles.py, RMS streamline
    # (60.8% span, interpolated between the 60%/70% tabulated sections), least-squares
    # fit of the SIGNED camber angle over the 65-88% chord window -- clear of the TE
    # FILLET that corrupted the last ~10% of chord and produced the wrong -32.95 deg
    # backsweep extraction (docs/centrifugal/11-blade-angle-discrepancy.md;
    # data/coefficients.md S1). Cross-checked by an independent 2-point secant over
    # the SAME window: -4.4779 rad/m (7% lower magnitude, same sign -- comparable
    # rigor to this project's other independent-read cross-checks, e.g. the TE metal
    # blockage's 0.3%/1% reads and l_blade_m's ~19% cross-check). Window sensitivity
    # (span-averaged, not RMS-span-specific): -4.01 to -4.82 rad/m across six window
    # choices from 55-80% to 75-88% chord, all clear of the fillet.
    #
    # THIS DISAGREES WITH THE PRE-REGISTERED FIGURE of -3.5 +- 1.5 rad/m
    # (docs/centrifugal/15-experiments-prereg.md E5) by about 30% in MAGNITUDE, though
    # it is inside the stated +-1.5 rad/m band (-5.0 to -2.0) and has the SAME SIGN --
    # the sign is what the falsifiable claim (Qiu S3.2; test_qiu_slip.py) depends on,
    # and it is not in question. Reported honestly rather than adjusted to the
    # pre-registered point value (docs/centrifugal/17-tdd-plan.md's own rule: "the
    # plan can be wrong").
    dbeta_dm_te: Optional[float] = None

    @property
    def Z_eff(self) -> float:
        """Splitter-aware effective blade count for slip (research/notes/
        01-centrifugal-meanline-methodology.md S3.6, S11; docs/centrifugal/
        02-tdd-plan.md).

        Aungier's simplified splitter treatment, as transcribed by Yang, Liu & Zhao
        2023 (*Machines* 11(1), 118), S1, attributed to Aungier [7] (citation
        corrected, docs/centrifugal/17-tdd-plan.md slice S0 -- the PDF's authors are
        Yang, Liu & Zhao, not "Li et al."; the file has since been renamed to
        ``research/papers/yang_2023_machines-11-118_loss-models-splitter.pdf`` to
        match, see research/papers/RENAME-MAP.md, verified):

            Z_eff = Z_FB + Z_SB * (L_SB / L_FB)

        with ``Z_FB``/``Z_SB`` the full/splitter blade counts and ``L_FB``/``L_SB``
        their MERIDIONAL LENGTHS (LE to TE, along the camberline).

        **Supply those lengths.** ``l_main_m``/``l_splitter_m`` take the formula's own
        quantity, and are then used verbatim. For HECC they are measured from NASA's
        blade coordinates (see the field comments).

        THE FALLBACK, AND WHY IT WAS WRONG
        ----------------------------------
        Without them, ``L`` is proxied by the RADIAL EXTENT the blade spans. The previous
        version of this proxy measured the full blade from the SHAFT CENTERLINE::

            L_FB = r_te            # from r = 0 -- where no blade exists
            L_SB = r_te - r_split

        A blade does not begin at the centerline; it begins at its leading edge. That
        inflates ``L_FB``, deflates ``L_SB/L_FB``, and biases ``Z_eff`` low. The proxy
        now spans LE to TE for BOTH blades, which is the same formula applied
        consistently::

            L_FB = r_te - r_le     # r_le from x_le, or the inducer-eye fallback
            L_SB = r_te - r_split

        Measured against NASA's actual blade camberlines (span-averaged, from
        ``my_scripts/extract_hecc_blade_angles.py``), the TRUE meridional ratio for HECC
        is ``L_SB/L_FB = 0.6937`` -> ``Z_eff = 25.405``. The old centerline proxy gave
        ``0.6872`` -> ``25.308``: conceptually indefensible but, as it happens, within 1%.
        So this correction is NOT what rescues (or breaks) any validation target -- it is
        fixed because it is wrong, not because it moves a number. Reported, either way.

        (An earlier audit estimated the true ``Z_eff`` at ~29.3 -- "13.5% low" -- from a
        flowpath arc-length estimate it flagged as unreliable in the near-axial inducer.
        Direct measurement of the blade camberlines refutes that: 25.4, not 29.3.)

        ``r_le`` for the proxy needs a leading edge. If ``x_le`` is not given, this falls
        back to ``splitter_le_r``-only information and finally to the naive
        ``Z_eff = Z_FB + Z_SB``. Callers who care must supply the lengths.
        """
        if self.n_splitters <= 0:
            return float(self.n_blades)

        # Preferred: the formula's own quantity, measured.
        if self.l_main_m is not None and self.l_splitter_m is not None:
            if self.l_main_m <= 0.0:
                raise ValueError("l_main_m must be > 0")
            return self.n_blades + self.n_splitters * (self.l_splitter_m / self.l_main_m)

        # Proxy: radial extent, LE -> TE for BOTH blades (never from r = 0).
        if self.splitter_le_r is None:
            return float(self.n_blades + self.n_splitters)
        r_le = 0.0 if self.r_le_hint is None else self.r_le_hint
        L_full = self.r_te - r_le
        L_split = self.r_te - self.splitter_le_r
        if L_full <= 0.0:
            return float(self.n_blades)
        return self.n_blades + self.n_splitters * (L_split / L_full)
