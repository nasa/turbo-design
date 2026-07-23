"""Slice 2 -- slip. Root cause C.

There is no slip factor anywhere in upstream turbodesign (grep for
slip|wiesner|busemann|stanitz returns only a vaneless-diffuser LOSS class that
happens to be named Stanitz). Without it the code assumes the blades guide the flow
perfectly, over-predicting Vt2 and hence work: psi = 0.955 against a measured 0.81.

Slip enters as a VELOCITY DEFICIT, which is its natural form:

    Vt2 = sigma*U2 + Cm2*tan(beta2b)          (beta2b < 0 for backsweep)

NOT as an additive deviation angle (upstream's DeviationBaseClass returns degrees,
which would force a lossy conversion).

WHAT THIS SLICE DOES *NOT* DO
-----------------------------
It does not make psi match NASA's measured 0.81. Slip alone gives the EULER work,
which lands at ~0.738. The remaining ~0.07 is PARASITIC work (disc friction,
recirculation, leakage) -- work that enters the fluid without producing pressure --
and it does not exist until Slice 4.

Do not tune sigma to close that gap. It is not slip's job to close it.
"""

import math

import pytest

from turbodesign.centrifugal import (
    Air,
    BusemannSlip,
    Impeller,
    InletState,
    MeridionalPath,
    PhysicsError,
    QiuSlip,
    StanitzSlip,
    Stage,
    WiesnerSlip,
)

R2 = 0.21581
RPM = 21789.0
MDOT = 4.9269
BETA2B = -37.5
N_BLADES = 15
N_SPLITTERS = 15


@pytest.fixture(scope="module")
def path(data_dir):
    return MeridionalPath.from_csv(
        data_dir / "hecc" / "flowpath_hub.csv",
        data_dir / "hecc" / "flowpath_shroud.csv",
    )


def _solve(path, slip, n_splitters=0):
    stage = Stage(
        path=path,
        impeller=Impeller(
            n_blades=N_BLADES,
            n_splitters=n_splitters,
            backsweep_deg=BETA2B,
            r_te=R2,
        ),
        diffuser=None,
        slip=slip,
        losses=None,
        fluid=Air(),
    )
    return stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))


# ------------------------------------------------- the behaviour change


def test_slip_reduces_the_work_factor(path):
    """The observable consequence. Slice 1 (no slip) gives psi = 0.872 -- perfect blade
    guidance, which is physically impossible. Slip brings it to ~0.738.
    """
    no_slip = _solve(path, slip=None)
    with_slip = _solve(path, slip=WiesnerSlip())

    assert no_slip.psi == pytest.approx(0.8719, rel=0.02)
    assert with_slip.psi == pytest.approx(0.7381, rel=0.02)
    assert with_slip.psi < no_slip.psi, "slip must REDUCE the work input"


def test_slip_does_not_close_the_gap_to_the_measured_work_factor(path):
    """DELIBERATE. NASA measures psi = 0.81; slip alone gives ~0.738.

    The residual ~0.07 is PARASITIC work and arrives in Slice 4. This test exists to
    stop anyone (including a future me) tuning sigma to hit 0.81 here -- that would
    silently bury the parasitic term inside the slip factor, and the two would then be
    inseparable forever.
    """
    op = _solve(path, slip=WiesnerSlip())
    assert op.psi < 0.79, "psi must still be BELOW the measured 0.81 -- parasitic work is missing"


# ------------------------------------------------- published slip factors (Tier C)
# Each model must reproduce the sigma its own source paper publishes.


def test_wiesner_slip_factor():
    """Wiesner (1967), Eq. 3:  sigma = 1 - sqrt(cos(beta2b)) / Z**0.7"""
    sigma = WiesnerSlip().sigma(beta2b_deg=BETA2B, Z=N_BLADES)
    expected = 1 - math.sqrt(math.cos(math.radians(abs(BETA2B)))) / N_BLADES**0.7
    assert sigma == pytest.approx(expected, rel=1e-9)
    assert sigma == pytest.approx(0.8662, abs=0.001)


def test_stanitz_slip_factor():
    """Stanitz (1952):  sigma = 1 - 0.63*pi/Z  (weak beta2b dependence)."""
    sigma = StanitzSlip().sigma(beta2b_deg=BETA2B, Z=N_BLADES)
    assert sigma == pytest.approx(0.868, abs=0.005)


def test_busemann_slip_factor_is_bounded():
    sigma = BusemannSlip().sigma(beta2b_deg=BETA2B, Z=N_BLADES)
    assert 0.80 < sigma < 0.95


@pytest.mark.parametrize("model", [WiesnerSlip(), StanitzSlip(), BusemannSlip(), QiuSlip()])
@pytest.mark.parametrize("Z", [8, 15, 20, 30])
@pytest.mark.parametrize("beta2b", [0.0, -20.0, -37.5, -50.0])
def test_slip_factor_is_always_physical(model, Z, beta2b):
    """Tier A: sigma is a velocity RATIO. It cannot exceed 1 or go negative."""
    sigma = model.sigma(beta2b_deg=beta2b, Z=Z, Cm2_over_U2=0.20)
    assert 0.0 < sigma <= 1.0


def test_more_blades_means_less_slip():
    """Monotonicity: better guidance -> sigma -> 1."""
    w = WiesnerSlip()
    sigmas = [w.sigma(beta2b_deg=BETA2B, Z=z) for z in (8, 15, 30, 60)]
    assert sigmas == sorted(sigmas), "sigma must increase with blade count"
    assert sigmas[-1] > 0.93


# ------------------------------------------------- Qiu: the default


def test_qiu_varies_with_operating_point():
    """Qiu et al. (2011) is the only model of the four here that is CAPABLE of
    responding to the operating point (via the exit flow coefficient, through the
    blade-turning term, Eq. 10b) -- which is why it remains the candidate for
    SPEEDLINES (slices 9-11), not just a design point. It is NOT the shipped default;
    Wiesner is (docs/centrifugal/17-tdd-plan.md S6; data/coefficients.md S5).

    Wiesner/Stanitz/Busemann are functions of geometry alone -- they return the same
    sigma at choke as at surge, which is wrong for a speedline.

    docs/centrifugal/17-tdd-plan.md slice S6: the phi2-dependence used to come from a
    wrong-signed PROXY that fired even with no ``dbeta_dm`` supplied. That proxy is
    DELETED -- absent ``dbeta_dm``, Qiu is now (correctly) geometry-only too, exactly
    like the other three (``tests/unit/test_qiu_slip.py::
    test_absent_gradient_makes_sigma_independent_of_flow_coefficient``). Supplying
    the real ``dbeta_dm``/``r2`` is what makes Qiu respond to phi2, per Eq. (10b).
    """
    q = QiuSlip()
    s_low = q.sigma(beta2b_deg=BETA2B, Z=N_BLADES, Cm2_over_U2=0.10, dbeta_dm=-3.5, r2=R2)
    s_high = q.sigma(beta2b_deg=BETA2B, Z=N_BLADES, Cm2_over_U2=0.35, dbeta_dm=-3.5, r2=R2)
    assert s_low != pytest.approx(s_high, rel=1e-6), "Qiu must respond to flow coefficient"

    w = WiesnerSlip()
    assert w.sigma(beta2b_deg=BETA2B, Z=N_BLADES, Cm2_over_U2=0.10) == pytest.approx(
        w.sigma(beta2b_deg=BETA2B, Z=N_BLADES, Cm2_over_U2=0.35), rel=1e-9
    ), "Wiesner is geometry-only, by construction"


# ------------------------------------------------- splitters


def test_splitters_reduce_slip(path):
    """HECC and CC3 are BOTH 15 main + 15 splitter. Ignoring splitters mis-predicts
    sigma materially: Z=15 -> 0.866, Z=30 -> 0.918.

    Z_eff must come from a CITED model, never from whichever value makes psi land on
    the measured 0.81. (Z_eff = 22.5 happens to hit it exactly. That is a trap, not a
    result -- see docs/centrifugal/04-acceptance-criteria.md.)
    """
    plain = _solve(path, slip=WiesnerSlip(), n_splitters=0)
    split = _solve(path, slip=WiesnerSlip(), n_splitters=N_SPLITTERS)
    assert split.psi > plain.psi, "splitters guide the flow better -> less slip -> more work"


def test_effective_blade_count_is_between_main_and_main_plus_splitter():
    """A splitter starts partway along the passage, so it cannot count as a full blade,
    and it certainly counts for more than nothing.
    """
    imp = Impeller(n_blades=15, n_splitters=15, backsweep_deg=BETA2B, r_te=R2, splitter_le_r=0.0675)
    assert 15 < imp.Z_eff < 30


# ------------------------------------------------- the guard that matters


def test_slip_on_a_vaned_diffuser_is_refused():
    """Qiu states this explicitly: with omega = 0 the Coriolis term vanishes and the
    slip derivation collapses. Use Carter's rule for a vaned diffuser.

    This is an easy, silent error to make once a SlipModel exists in the codebase --
    so it must RAISE, not quietly compute a meaningless number.
    """
    with pytest.raises(PhysicsError, match="(?i)coriolis|rotat|slip"):
        QiuSlip().sigma(beta2b_deg=-10.0, Z=20, omega=0.0)
