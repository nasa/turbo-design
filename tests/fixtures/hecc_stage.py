"""HECC validation driver — the whole stage, on the new centrifugal module.

    uv run python tests/fixtures/hecc_stage.py

A test fixture, not library API: it builds NASA's HECC stage and solves it, and is
imported by tests/unit and tests/oracle for the HECC design point. Also runnable
standalone (see above). Replaces the old nasa-hecc.py, which drove the UPSTREAM
axial-oriented solver and therefore carried every defect that solver has (see
docs/centrifugal/00-root-cause-analysis.md). This one drives turbodesign.centrifugal
end to end: inducer -> impeller -> vaneless space -> vaned diffuser -> exit guide vanes.

EVERY GEOMETRIC NUMBER BELOW IS SOURCED. Nothing is fitted, and nothing is tuned to make
the answer look better. Where the model misses, it misses, and the miss is printed.

Outputs
-------
  tests/fixtures/outputs/hecc_speedline.png   PR and eta vs mass flow, model vs NASA measured
  tests/fixtures/outputs/hecc_summary.txt     the design-point table and the honest deltas
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

from turbodesign.centrifugal import (
    Air,
    Impeller,
    InletState,
    MeridionalPath,
    OhLossSet,
    Stage,
    WiesnerSlip,
)
from turbodesign.centrifugal.diffusion import (
    ExitGuideVanes,
    VanedDiffuser,
    VanelessSpace,
)

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "outputs"
IN = 0.0254  # in -> m

# ---------------------------------------------------------------- the machine, sourced
#
# NASA/CR-2014-218114/REV1 (Medic et al.), Table 2, unless noted.
RPM = 21789.0  # data/hecc/design_point.csv
MDOT_DESIGN = 4.9269  # kg/s corrected, ibid.
P01, T01 = 101325.0, 288.15  # standard day; mdot is the CORRECTED flow

# NASA's measured design point (data/hecc/design_point.csv) -- the targets.
NASA = {"psi": 0.81, "PR_tt": 4.6847, "eta_poly": 0.8553}

IMPELLER = dict(
    n_blades=15,
    n_splitters=15,
    r_te=0.21581,  # 16.988/16.998 in exit dia
    splitter_le_r=0.0675,
    tip_clearance=0.000305,  # 0.012 in. HECC's measured schedule is FLAT
    # (0.010/0.0125/0.012 in from 3 tip-capacitance probes) -- the 4x schedule in the
    # literature belongs to CC3, not HECC. See docs/centrifugal/09-tip-clearance.md.
    # --- derived from NASA's own Appendix C blade coordinates, NOT tabulated by NASA ---
    inducer_blade_angle_deg=45.46,  # RMS streamline (60.8% span), blade_angles.csv
    l_main_m=0.209875,  # true meridional camberline length
    l_splitter_m=0.145580,
    # TRUE through-blade (camberline) length -- docs/centrifugal/17-tdd-plan.md slice
    # S1. Feeds skin friction + mixing (NOT leakage, which uses l_main_m above -- a
    # DIFFERENT quantity).
    #
    # MEASURED, not estimated: the 3-D camberline arc ds = sqrt(dm^2 + (r*dtheta)^2),
    # integrated along each of NASA Appendix C's 11 main-blade sections and span-averaged
    # (extract_hecc_blade_angles.py). It needs NO blade angle at all -- it is
    # the geometry.
    #
    # NOT l_main_m / cos(beta1b) = 0.292395 m (+23%). That divides the WHOLE blade's
    # meridional length by the cosine of the INDUCER LEADING-EDGE angle (45.46 deg) --
    # a one-station quantity applied to the whole blade. HECC's blade is S-shaped, so
    # beta is nowhere near its LE value over most of the chord (NASA Fig 18: a mid-chord
    # minimum near 14-20 deg, cos ~ 0.95, not 0.70). Sixth instance of this project's
    # signature failure mode -- see docs/PHYSICS-RULES.md.
    #
    # Replaces the shipped closure's L_tilde = 0.522 m, which used the LE flow
    # coefficient (Vm1/U1 = 0.759) where the source's own third term is the GLOBAL
    # flow coefficient (0.0608) -- 12.5x too large in the dominant term (data/coefficients.md).
    l_blade_m=0.237875,
    throat_area=0.020525,  # min-distance between adjacent MAIN blades; 15 passages
    # LE TANGENTIAL BLADE THICKNESS, for Conrad's blockage-corrected optimum-incidence
    # angle. Read the definition before you trust the number.
    #
    # NASA's tabulated blade surfaces MEET AT A POINT at the LE: t(m=0) is EXACTLY 0, and
    # the blade then thickens like a wedge with NO PLATEAU (t = 0.100 in at 0.5% chord,
    # 0.159 at 5%, 0.169 at 10%). So "the LE thickness" is not a measurement -- it is a
    # CHOICE OF STATION, and it is a free knob of exactly the kind this project has been
    # burned by: across that range it moves the incidence loss 65x (260 -> 4 J/kg) while
    # moving psi by 0.001 and eta by 0.0006. NO TEST CAN SEE IT.
    #
    # It is pinned to the INDUCER THROAT -- the station where the passage blockage the
    # incoming flow accelerates through is actually established, and the one near-LE
    # station we have ALREADY located from NASA's coordinates (it is the same cut that
    # gives throat_area above). At the RMS streamline the throat sits at 7.98% of
    # meridional chord. extract_hecc_le_thickness.py.
    le_blade_thickness=0.00418982,  # 0.1650 in; blockage Z*t/(pi*D1) = 12.3%
    # EXIT (TE) METAL BLOCKAGE -- docs/centrifugal/17-tdd-plan.md slice S4. SOLID METAL,
    # not the aerodynamic (boundary-layer/wake) blockage this project already killed
    # once (data/coefficients.md) -- see Impeller.blockage's own docstring.
    #
    # Both blade rows (15 main + 15 splitter) reach the trailing edge and both have a
    # nonzero TANGENTIAL thickness there -- exactly like the vaned diffuser's own
    # already-accepted TE-thickness blockage below (blockage=0.031). Derived from NASA
    # Appendix C by the SAME method already accepted for the LE thickness t1
    # (extract_hecc_le_thickness.py):
    #
    #   B_geom = (Z_main*t_main_TE + Z_splitter*t_splitter_TE) / (2*pi*r2)
    #          = (15*0.032200 + 15*0.024164) / (2*pi*8.4965)   [inches]
    #          = 0.015837
    #
    # STATION CHOICE, exactly like t1's: NASA's tabulated surfaces close to a
    # near-point at the TE (t_theta ~ 0.032 in at 100% chord, but 0.138 in at 99% and
    # 0.171 in at 98% -- a >5x lever). 100% chord (the true exit plane) is ADOPTED
    # because there is no metal downstream of r2 -- not because of what it does to the
    # answer. extract_hecc_te_blockage.py; independent (quadratic-
    # interpolation) cross-check agrees to 0.3%.
    blockage=0.015837,
    # SIGNED (dbeta/dm)_2 at the TE, rad/m -- docs/centrifugal/17-tdd-plan.md slice
    # S6. Feeds QiuSlip's dsigma_turn (Qiu 2011 Eq. 10b) ONLY -- WiesnerSlip (the
    # slip model configured below) ignores this field entirely via its **kw seam, so
    # populating it here is OUTPUT-INVARIANT (tests/unit/test_qiu_slip.py::
    # test_shipped_outputs_are_bit_identical).
    #
    # DERIVED, not the plan's pre-registered -3.5 rad/m: least-squares fit of the
    # signed camber angle over the fillet-clear 65-88% chord window, at the RMS
    # streamline (60.8% span), extract_hecc_blade_angles.py. Cross-checked
    # by an independent 2-point secant over the same window: -4.4779 rad/m (7% lower
    # magnitude, same sign). See Impeller.dbeta_dm_te's own docstring for the window
    # sensitivity and the honest disagreement with the pre-registered figure.
    dbeta_dm_te=-4.8255,
)

# THE EXIT BLADE ANGLE -- settled, and NOT by us choosing.
#
# History, because it is the whole methodological point of this case. Three values have
# stood here:
#   -37.5  UNSOURCED. The arithmetic midpoint of NASA's PROSE band ("32 to 42 deg"). It was
#          the ONLY value at which the psi gate passed -- which is how it survived.
#   -32.95 DERIVED, and WRONG. Span-average of the TE metal angle extracted from NASA's own
#          Appendix C coordinates. Correct machinery, corrupted input: the impeller has a
#          ROUNDED trailing edge (Fig 17), Appendix C tabulates the surface loop, so the
#          loop WRAPS THE FILLET and the camber angle rolls over in the last ~10% of chord.
#          The extraction's TE fit window (85-95% chord) sits INSIDE that rolled-over zone.
#   -36.0  MEASURED BY NASA. Figure 18 of CR-2014-218114/REV1 (PDF p. 36) plots exactly this
#          quantity -- blade angle vs %chord at 5 spans -- for this impeller. Digitized at
#          400 dpi, +-0.5 deg: TE = 31.3/32.9/35.2/38.2/42.5 deg at 0/25/50/75/100% span.
#          At a RADIAL TE, dA = 2*pi*r*dx with r constant, so the area-average IS the plain
#          span-average: 36.0 deg. See docs/centrifugal/11-blade-angle-discrepancy.md.
#
# Fig 18 also CONFIRMS our leading-edge angles to ~1 deg, which refutes the alternative
# hypothesis (that NASA used a different angle convention) -- a convention change would have
# moved the LE too. The error was localised to the fillet, and it is now gone.
#
# psi has slope 0.0066/deg here, so this is the most leveraged number in the model. It is
# now the one number in it that NASA measured and published a plot of.
BACKSWEEP = -36.0  # NASA CR-2014-218114/REV1 Fig 18, span-averaged. +-0.5 deg read, ~1 deg span-weighting.


def diffusion_system() -> list:
    """Vaneless space -> vaned diffuser -> EGV. Table 2 + Appendix C vane angles."""
    return [
        VanelessSpace(r3=9.1074 * IN, b3=0.559 * IN),
        VanedDiffuser(
            r3=9.1074 * IN,  # 18.215 in dia
            r4=11.2044 * IN,  # 22.398 in dia
            b=0.559 * IN,  # channel height (a PINCH from the impeller's 0.609)
            n_vanes=20,
            n_splitters=20,
            beta_le_deg=79.64,  # Appendix C C.25/C.26
            beta_te_deg=36.71,
            chord=2.097 * IN,
            # SPLITTER VANE chord + LE radius -- docs/centrifugal/17-tdd-plan.md slice S5.
            # The splitter LE sits at r=0.2488 m, INSIDE the r3=0.23133 -> r4=0.28459 m
            # span -- only 20 (not 40) vanes exist over the first third of the passage,
            # exactly where circulation per vane is largest. These feed ONLY
            # VanedDiffuser.Z_eff_loading (the diffusion factor / loading loss); the
            # deviation and friction cascade above keeps ALL 40 vanes (a TE / wetted-
            # passage quantity -- both blade rows reach the TE).
            #
            # MEASURED, not estimated (extract_hecc_vane_angles.py, from NASA
            # Appendix C Tables C.27/C.28): chord_splitter = 1.4110 in = 0.0358396 m,
            # splitter_le_r = 9.7952 in = 0.248798 m -> Z_eff,vd = 33.4575.
            chord_splitter=0.0358396286,
            splitter_le_r=0.248797826,
            blockage=0.031,  # GEOMETRIC vane metal blockage, from the TE thicknesses
            throat_area=0.005983,  # 20 passages; splitters do NOT block it
        ),
        ExitGuideVanes(
            r_mean=12.056 * IN,
            b=0.4803 * IN,
            n_vanes=60,
            chord=2.434 * IN,
            beta_le_deg=48.1,  # Appendix C C.29/C.30
            beta_te_deg=3.2,
        ),
    ]


def build(backsweep: float) -> Stage:
    path = MeridionalPath.from_csv(
        REPO / "data" / "hecc" / "flowpath_hub.csv",
        REPO / "data" / "hecc" / "flowpath_shroud.csv",
    )
    return Stage(
        path=path,
        impeller=Impeller(backsweep_deg=backsweep, **IMPELLER),
        slip=WiesnerSlip(),
        losses=OhLossSet(),
        fluid=Air(),
        components=diffusion_system(),
    )


def nasa_speedline() -> list[tuple[float, float, float]]:
    with open(REPO / "data" / "hecc" / "map_vaned.csv") as fh:
        rows = list(csv.DictReader(fh))
    pts = [
        (float(r["mdot_corr_kgs"]), float(r["PR_tt"]), float(r["eta_poly"]))
        for r in rows
    ]
    return sorted(pts)


def design_point_table() -> str:
    """The design point, on the one geometry NASA published a plot of.

    eta is reported REAL-GAS. The constant-cp polytropic exponent inflates eta by ~0.8
    points, because eta is DEFINED through gamma and our internal solve freezes the gas at
    288 K. A documented approximation is still an error when the measurement you compare
    against did not make it.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        op = build(BACKSWEEP).solve(
            mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )

    lines = [
        "HECC DESIGN POINT -- model vs NASA measured",
        "",
        f"backsweep = {BACKSWEEP} deg (NASA Fig 18, digitized). Every geometric input is",
        "sourced. No coefficient is tuned. Where it misses, it misses, and the miss is printed.",
        "",
        f"{'':24s}{'model':>10}{'NASA':>10}{'error':>10}",
        "-" * 54,
    ]
    got = {
        "psi": op.psi,
        "PR_tt": op.stage_PR,
        "eta_poly": op.stage_eta_poly_realgas,
    }
    for key, name in (
        ("psi", "work factor psi"),
        ("PR_tt", "stage PR_tt"),
        ("eta_poly", "stage eta_poly"),
    ):
        a, n = got[key], NASA[key]
        lines.append(f"{name:24s}{a:>10.4f}{n:>10.4f}{100 * (a / n - 1):>+9.1f}%")

    lines += [
        "",
        f"(constant-cp eta_poly would read {op.stage_eta_poly:.4f} -- inflated by "
        f"{100 * (op.stage_eta_poly - op.stage_eta_poly_realgas):.1f} points. Not reported.)",
        "",
        "WHAT IS STILL MISSING -- read before quoting any number above.",
        "",
        "  1. The vaned diffuser runs at a diffusion factor of Df ~= 1.14, well beyond the",
        "     fitted range of EVERY published cascade correlation (Lieblein's own validity",
        "     limit is Df < 0.6). An earlier version of this line additionally claimed that",
        "     above Df >= 1 'Lieblein's wake-thickness law contains ln(1 - Df) and is",
        "     UNDEFINED' -- that specific functional form is RETRACTED as UNVERIFIED (slice",
        "     S0, docs/centrifugal/17-tdd-plan.md): NACA RM E53D01 is not held and the claim",
        "     could not be checked against the primary. What stands regardless: every",
        "     correlation descended from Lieblein stops being fitted near Df ~ 0.6, so there",
        "     is no published loss law at Df ~= 1.14 either way, and we have no sourced",
        "     blockage model for a separated diffuser. That is worth roughly 1 eta point, and",
        "     it is NOT filled -- inventing a blockage coefficient here would be exactly the",
        "     free knob this project has been burned by four times.",
        "",
        "  2. CAPACITY, not loss: the model chokes ~11% high in mass flow (see the speedline",
        "     -- model chokes near 5.82 kg/s, NASA chokes at 5.24). No loss coefficient",
        "     appears in the choke-flow expression -- it is set by throat area, P02/sqrt(T02)",
        "     and gamma alone -- so this is a separate, loss-independent error and cannot be",
        "     absorbed by any tuning.",
        "",
        "  3. The internal solve still uses constant cp (worth ~1.5 K on T02). The eta above",
        "     is real-gas corrected at the reporting boundary, not in the solve.",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    table = design_point_table()
    print(table)
    (OUT / "hecc_summary.txt").write_text(table + "\n")

    # ---- speedline
    nasa = nasa_speedline()
    stage = build(BACKSWEEP)
    flows, prs, etas = [], [], []
    m = 4.60
    while m <= 5.90:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                op = stage.solve(mdot=m, rpm=RPM, inlet=InletState(P0=P01, T0=T01))
            if op.stage_PR and op.stage_PR > 1.0:
                flows.append(m)
                prs.append(op.stage_PR)
                etas.append(op.stage_eta_poly_realgas)  # NOT the constant-cp value
        except Exception:  # choked / no subsonic solution -- the map simply ends
            break
        m += 0.02

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib not installed -- skipping the plot)")
        return

    NASA_CHOKE = 5.24
    m_choke = flows[-1]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dp = stage.solve(mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))
    nm = [p[0] for p in nasa]

    for ax, model_y, nasa_y, dp_y, ylab, title in (
        (
            ax1,
            prs,
            [p[1] for p in nasa],
            dp.stage_PR,
            "stage $PR_{tt}$",
            "Pressure ratio",
        ),
        (
            ax2,
            etas,
            [p[2] for p in nasa],
            dp.stage_eta_poly_realgas,
            r"stage $\eta_{poly}$",
            "Polytropic efficiency (real-gas)",
        ),
    ):
        ax.axvline(NASA_CHOKE, color="0.55", ls=":", lw=1.4, zorder=1)
        ax.axvline(m_choke, color="crimson", ls=":", lw=1.4, zorder=1)
        ax.plot(nm, nasa_y, "ko-", ms=5, lw=1.6, label="NASA measured", zorder=3)
        ax.plot(
            flows,
            model_y,
            "-",
            color="crimson",
            lw=2.2,
            label="meanline model",
            zorder=3,
        )
        ax.plot(
            [MDOT_DESIGN],
            [dp_y],
            "*",
            color="crimson",
            ms=15,
            mec="k",
            mew=0.6,
            label="model, design point",
            zorder=4,
        )
        ax.set_xlabel("corrected mass flow, kg/s")
        ax.set_ylabel(ylab)
        ax.set_title(title, fontsize=11)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower left")

    # the capacity error, annotated on the PR panel -- it is the headline miss
    y_arrow = 3.15
    ax1.annotate(
        "",
        xy=(m_choke, y_arrow),
        xytext=(NASA_CHOKE, y_arrow),
        arrowprops=dict(arrowstyle="<->", color="crimson", lw=1.4),
    )
    ax1.text(
        0.5 * (NASA_CHOKE + m_choke),
        y_arrow + 0.10,
        f"+{100 * (m_choke / NASA_CHOKE - 1):.0f}% capacity error\n(loss-INDEPENDENT)",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="crimson",
    )
    ax1.text(
        NASA_CHOKE - 0.04,
        4.62,
        "NASA\nchokes\n5.24",
        ha="right",
        va="center",
        fontsize=8,
        color="0.35",
    )

    fig.suptitle(
        f"NASA HECC, 100% speed — meanline model vs measurement.   "
        f"backsweep {BACKSWEEP}° (NASA Fig 18); no coefficient tuned.\n"
        f"Design point: ψ {dp.psi:.4f} (+{100 * (dp.psi / NASA['psi'] - 1):.1f}%),  "
        f"PR {dp.stage_PR:.3f} (+{100 * (dp.stage_PR / NASA['PR_tt'] - 1):.1f}%),  "
        f"η {dp.stage_eta_poly_realgas:.4f} (+{100 * (dp.stage_eta_poly_realgas / NASA['eta_poly'] - 1):.1f}%)"
        f"   —   the SHAPE is right, the CAPACITY is not.",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "hecc_speedline.png", dpi=150)
    print(f"\nwrote {(OUT / 'hecc_speedline.png').relative_to(REPO)}")
    print(f"      model chokes near {m_choke:.2f} kg/s; NASA chokes at {NASA_CHOKE}")


if __name__ == "__main__":
    main()
