"""Frozen independent-calculation oracle for NASA HECC.

This does NOT use turbodesign. It is a from-scratch 1D meanline calculation with
correct physics -- rothalpy-consistent P0R and Wiesner slip -- built to answer one
question: does the physics close at all?

It does. With UN-TUNED, guessed loss coefficients it reaches:

    psi = 0.751   PR = 4.13   eta_poly = 83.1%     (NASA: 0.80 / 4.68 / 85.5%)

while turbodesign, on identical geometry and the same design point, returns:

    psi = 0.955   PR = 2.99   eta_poly = 55.1%

That gap is the whole bug. See docs/centrifugal/00-root-cause-analysis.md.

These numbers are FROZEN. `turbodesign.centrifugal` must at minimum match this oracle
(slice 1 asserts the isentropic impeller PR of ~5.35 that falls out of it), and with a
calibrated Oh-1997 loss set should beat it.

Do not "improve" the oracle to make the library agree with it. If they disagree, one of
them is wrong and the disagreement is the finding.
"""

import numpy as np
import pytest

CP, GAMMA, R_GAS = 1005.0, 1.4, 287.0

# --- NASA HECC design point (NASA/CR-2014-218114/REV1, Table 2) ---
MDOT, RPM = 4.95, 21789.0
T01, P01 = 288.15, 101325.0
R1_HUB, R1_SHROUD = 0.04051, 0.10770
R2, B2 = 0.21581, 0.01547
N_BLADES, BETA2B_DEG = 15, -37.5  # backsweep, negative = backswept


def _solve_inlet(cx_lo=50.0, cx_hi=250.0):
    """Axial inlet velocity from continuity (bisection)."""
    area = np.pi * (R1_SHROUD**2 - R1_HUB**2)

    def residual(cx):
        t1 = T01 - cx**2 / (2 * CP)
        p1 = P01 * (t1 / T01) ** (GAMMA / (GAMMA - 1))
        return p1 / (R_GAS * t1) * cx * area - MDOT

    lo, hi = cx_lo, cx_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(lo) * residual(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def hecc_oracle() -> dict:
    omega = RPM * np.pi / 30
    u2 = omega * R2

    # STREAMLINE CONSISTENCY: the relative-frame state must be evaluated on ONE
    # streamline. We use the RMS radius (the standard meanline choice -- it splits
    # the inlet annulus into equal areas). The shroud is carried separately, purely
    # as an inducer-Mach diagnostic; it must never enter the relative-frame closure.
    #
    # (An earlier version of this oracle built T0R1 from the SHROUD relative velocity
    # while taking the exit state at the meanline. Mixing streamlines like that gave
    # P0R2,is/P0R1 = 2.40 instead of 2.78. The unit test caught it.)
    r1_rms = np.sqrt(0.5 * (R1_HUB**2 + R1_SHROUD**2))
    u1 = omega * r1_rms
    u1_shroud = omega * R1_SHROUD

    # Wiesner slip
    sigma = 1 - np.sqrt(np.cos(np.radians(abs(BETA2B_DEG)))) / N_BLADES**0.7

    cx1 = _solve_inlet()
    t1 = T01 - cx1**2 / (2 * CP)
    p1 = P01 * (t1 / T01) ** (GAMMA / (GAMMA - 1))

    w1 = np.hypot(cx1, u1)  # meanline relative velocity -- used for the closure
    w1s = np.hypot(cx1, u1_shroud)  # shroud -- DIAGNOSTIC ONLY
    m1s_rel = w1s / np.sqrt(GAMMA * R_GAS * t1)

    t0r1 = t1 + w1**2 / (2 * CP)
    p0r1 = p1 * (t0r1 / t1) ** (GAMMA / (GAMMA - 1))

    # impeller exit: iterate Cm2 on continuity
    cm2 = 0.25 * u2
    for _ in range(300):
        vt2 = sigma * u2 + cm2 * np.tan(np.radians(BETA2B_DEG))
        work = u2 * vt2  # Euler, no inlet swirl
        t02 = T01 + work / CP
        w2 = np.hypot(cm2, vt2 - u2)
        t2 = t02 - (cm2**2 + vt2**2) / (2 * CP)
        t0r2 = t2 + w2**2 / (2 * CP)

        # THE term turbodesign discards (compressor_math.py:246):
        p0r2_is = p0r1 * (t0r2 / t0r1) ** (GAMMA / (GAMMA - 1))

        dh_incidence = 0.5 * (w1 * np.sin(np.radians(12))) ** 2
        dh_skin = 2 * 0.005 * (0.5 * (w1 + w2)) ** 2 * 2.0
        dh_loading = 0.05 * work * 0.25
        dh_clearance = (
            0.6 * (0.0003 / B2) * abs(vt2) * np.sqrt(abs(4 * np.pi / (B2 * N_BLADES) * (vt2 * cm2)))
        )
        dh_internal = dh_incidence + dh_skin + dh_loading + dh_clearance

        p0r2 = p0r2_is * np.exp(-dh_internal / (R_GAS * t2 * GAMMA / (GAMMA - 1)))
        p2 = p0r2 / (t0r2 / t2) ** (GAMMA / (GAMMA - 1))
        rho2 = p2 / (R_GAS * t2)
        cm2_new = MDOT / (rho2 * 2 * np.pi * R2 * B2 * 0.95)
        if abs(cm2_new - cm2) < 1e-8:
            break
        cm2 = 0.5 * cm2 + 0.5 * cm2_new

    p02 = p2 * (t02 / t2) ** (GAMMA / (GAMMA - 1))

    # parasitic: adds work, no pressure
    dh_parasitic = (0.02 + 0.015 + 0.01) * work
    t02_actual = T01 + (work + dh_parasitic) / CP
    psi = CP * (t02_actual - T01) / u2**2

    # diffuser + deswirler: static recovery, then residual dynamic head
    c2 = np.hypot(cm2, vt2)
    q2 = 0.5 * rho2 * c2**2
    p4 = p2 + 0.60 * q2
    p04 = min(p4 + 0.5 * q2 * 0.18, p02 * 0.93)

    pr = p04 / P01
    t04 = t02_actual
    eta_is = (pr ** ((GAMMA - 1) / GAMMA) - 1) / (t04 / T01 - 1)
    eta_poly = (GAMMA - 1) / GAMMA * np.log(pr) / np.log(t04 / T01)

    return {
        "U2": u2,
        "sigma": sigma,
        "Cm2": cm2,
        "Vt2_over_U2": vt2 / u2,
        "psi": psi,
        "PR_impeller": p02 / P01,
        "PR": pr,
        "eta_is": eta_is,
        "eta_poly": eta_poly,
        "P0R2_is_over_P0R1": p0r2_is / p0r1,
        "M1_shroud_rel": m1s_rel,
    }


@pytest.fixture(scope="module")
def oracle():
    return hecc_oracle()


def test_tip_speed(oracle):
    assert oracle["U2"] == pytest.approx(492.0, abs=1.0)


def test_wiesner_slip_factor(oracle):
    assert oracle["sigma"] == pytest.approx(0.866, abs=0.005)


def test_the_discarded_pressure_term_is_large(oracle):
    """turbodesign assumes this ratio is exactly 1.0 (compressor_math.py:246).

    It is not. Discarding it is the entire 54% efficiency hole.

    The exact value depends on which inlet streamline the relative frame is
    evaluated on -- 2.87 at the arithmetic midspan (what turbodesign's single
    streamline uses), 2.78 at the RMS radius (what this oracle uses), 2.40 at the
    shroud. The LOAD-BEARING claim is not the third digit; it is that the factor is
    >> 1.0 on any streamline.
    """
    assert oracle["P0R2_is_over_P0R1"] == pytest.approx(2.78, abs=0.05)
    assert oracle["P0R2_is_over_P0R1"] > 2.0, "if this is ~1.0, the diagnosis is wrong"


def test_work_factor_is_physical(oracle):
    """turbodesign returns 0.955 (no slip). NASA measures 0.79-0.81."""
    assert oracle["psi"] == pytest.approx(0.751, abs=0.02)
    assert oracle["psi"] < 0.85, "psi ~0.95 means slip is missing"


def test_stage_pressure_ratio(oracle):
    """Un-tuned losses. NASA measures 4.68; turbodesign returns 2.99."""
    assert oracle["PR"] == pytest.approx(4.13, abs=0.10)
    assert oracle["PR"] > 3.5, "must beat turbodesign's 2.99 by a wide margin"


def test_polytropic_efficiency(oracle):
    """NASA measures 85.5%; turbodesign returns 55.1%."""
    assert oracle["eta_poly"] == pytest.approx(0.831, abs=0.01)
    assert oracle["eta_poly"] > 0.75, "eta ~0.55 means work is being destroyed"


def test_exit_flow_coefficient_is_physical(oracle):
    """Cm2/U2 ~ 0.20-0.30 for a real centrifugal.

    With b2 = 0 (the original geometry bug) this would need Cm2 = 11,042 m/s.
    """
    assert 0.15 < oracle["Cm2"] / oracle["U2"] < 0.30
