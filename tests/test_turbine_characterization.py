# tests/test_turbine_characterization.py
"""Pin the current behaviour of the axial and radial turbine examples.

Why this file exists: `flow_math.compute_streamline_areas` computes streamtube areas with
`2*pi*C*(S/2*dx**2 + r*dx)`, where S is a radius difference in metres. The first term is a
cubic metre and the second a square metre, so the expression adds a volume to an area (see
tests/test_band_area_oracle.py). That function is shared -- `turbine_math.py:306` calls it
directly and `turbine_spool.py` reaches it through `compute_massflow` -- so correcting the
area formula will change axial turbine results too, and axial turbines are this library's
primary published use case.

These are characterization goldens: they record what the code does today, not what is
physically correct. `eee_hpt.py` is expected to move once the shared area fix lands; the
point of pinning it now is that the change is then visible rather than discovered.
`radial_turbine-1D.py` is pinned because it should move further: it is a radial machine, so
it takes the buggy branch on every band. Over its solve, `compute_streamline_areas` is
called 64 times and the radial branch body executes 256 times (64 calls x 4 bands at 5
streamlines), at meridional angles up to 89.96 degrees. The axial branch never runs.

Both sets of goldens were produced by running the examples on this tree (see
`tests/conftest.py::_run`).
"""

import pytest

RTOL = 1e-9


# ---------------------------------------------------------------------------
# examples/EEE-HPT/eee_hpt.py -- 2-stage axial turbine, SLSQP multi-stage pressure
# balance (turbodesign/turbine_spool.py:_balance_pressure, fmin_slsqp branch).
#
# This example has no main() -- it is a flat top-level script -- so `example_spool`
# returns its module namespace directly. It defines `stator1`, `rotor1`, `stator2`,
# `rotor2` (the same BladeRow objects that are inside `spool.rows`, mutated in place
# by `spool.solve()`) and `spool` itself. Runtime is ~7-8s (SLSQP + 5 streamlines),
# hence @pytest.mark.slow.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_eee_hpt_is_five_streamlines(example_spool):
    spool = example_spool("EEE-HPT/eee_hpt.py")["spool"]
    assert spool.num_streamlines == 5


@pytest.mark.slow
def test_eee_hpt_massflow_attribute_is_the_input_echo(example_spool):
    """`spool.massflow` is set once at construction (turbine_spool.py:92) and is never
    written back to during `.solve()` in pressure-balance mode; only
    `solve_massflow_for_power`, which this example does not use, updates it. So the
    attribute is not the converged massflow -- it is the 20 kg/s guess the example passed
    in. The converged value is pinned separately below, from `spool.convergence_history`.
    """
    spool = example_spool("EEE-HPT/eee_hpt.py")["spool"]
    assert spool.massflow == pytest.approx(20.0, rel=RTOL)


@pytest.mark.slow
def test_eee_hpt_stator1_per_streamline(example_spool):
    row = example_spool("EEE-HPT/eee_hpt.py")["stator1"]
    assert row.P0 == pytest.approx([1230607.077101921] * 5, rel=RTOL)
    assert row.T0 == pytest.approx([1587.0] * 5, rel=RTOL)
    assert row.M == pytest.approx(
        [
            0.9008453280523004,
            0.8753927855297182,
            0.8517684601550825,
            0.8297668414327242,
            0.809148370680553,
        ],
        rel=RTOL,
    )
    assert row.Yp == pytest.approx([0.057] * 5, rel=RTOL)


@pytest.mark.slow
def test_eee_hpt_rotor1_per_streamline(example_spool):
    row = example_spool("EEE-HPT/eee_hpt.py")["rotor1"]
    assert row.P0 == pytest.approx(
        [
            578345.3532063541,
            568953.3631594869,
            561067.2952288443,
            554464.6265053176,
            549027.9543868615,
        ],
        rel=RTOL,
    )
    assert row.T0 == pytest.approx(
        [
            1343.4794419417044,
            1338.351413116741,
            1334.0491926457082,
            1330.470258091159,
            1327.5634278022285,
        ],
        rel=RTOL,
    )
    assert row.M == pytest.approx(
        [
            0.39115615492041317,
            0.3761548794976443,
            0.36377563900765325,
            0.35378845466081904,
            0.3458730120679529,
        ],
        rel=RTOL,
    )
    assert row.Yp == pytest.approx([0.088] * 5, rel=RTOL)


@pytest.mark.slow
def test_eee_hpt_stator2_per_streamline(example_spool):
    row = example_spool("EEE-HPT/eee_hpt.py")["stator2"]
    assert row.P0 == pytest.approx(
        [
            564286.1117135661,
            555542.1689799328,
            548200.2397365044,
            542053.1551549011,
            536991.6134126185,
        ],
        rel=RTOL,
    )
    assert row.M == pytest.approx(
        [
            0.8623745634810988,
            0.8197532030358293,
            0.7825474243432641,
            0.7498246933734738,
            0.7206268276705512,
        ],
        rel=RTOL,
    )
    assert row.Yp == pytest.approx([0.069] * 5, rel=RTOL)


@pytest.mark.slow
def test_eee_hpt_rotor2_per_streamline(example_spool):
    row = example_spool("EEE-HPT/eee_hpt.py")["rotor2"]
    assert row.P0 == pytest.approx(
        [
            270414.06537260045,
            263897.33914370113,
            259578.78860271804,
            257046.7352292148,
            256108.09563705017,
        ],
        rel=RTOL,
    )
    assert row.T0 == pytest.approx(
        [
            1134.8998922573246,
            1128.4000244077208,
            1124.1007855487087,
            1121.6582500344703,
            1120.9101696617186,
        ],
        rel=RTOL,
    )
    assert row.M == pytest.approx(
        [
            0.4721947561743981,
            0.4490594317991563,
            0.43211060563523745,
            0.420354573458027,
            0.4127988942866575,
        ],
        rel=RTOL,
    )
    assert row.Yp == pytest.approx([0.014] * 5, rel=RTOL)


@pytest.mark.slow
def test_eee_hpt_slsqp_tripwire(example_spool):
    """A coarse tripwire, deliberately looser than the row-level goldens above.

    `fmin_slsqp`'s own `x`/`fun` live inside `TurbineSpool._balance_pressure`
    (turbodesign/turbine_spool.py) and are never assigned to an attribute the example can
    see, so they cannot be pinned without reaching into library internals.
    `spool.convergence_history[-1]` records the same optimizer's final objective value
    (`massflow_std`, what `fmin_slsqp` was minimizing) and its iteration count, which
    answers the same question: did the solver converge to roughly the same place. For
    anything that needs to be exact, prefer the per-row assertions above.
    """
    spool = example_spool("EEE-HPT/eee_hpt.py")["spool"]
    assert len(spool.convergence_history) == 106
    last = spool.convergence_history[-1]
    assert last["massflow_std"] == pytest.approx(0.0014191633426164, rel=1e-6)
    assert last["massflow"] == pytest.approx(29.2606442034497, rel=1e-6)


# ---------------------------------------------------------------------------
# examples/radial-turbine/radial_turbine-1D.py -- single-stage radial-inflow turbine,
# single-stage pressure balance (turbodesign/turbine_spool.py:_balance_pressure,
# minimize_scalar branch -- only one stage, so fmin_slsqp is not used here).
#
# Pinned because it is expected to move once the area fix lands: this passage takes the
# `S/2 * dx**2` branch of `compute_streamline_areas` (flow_math.py:40) on every band --
# 256 band evaluations over 64 calls -- at meridional angles up to 89.96 degrees, where
# the term being added is dimensionally a volume rather than an area. On an axial machine
# that term is small; here it is not. This example also has no main(); `example_spool`
# returns its module namespace, which defines `stator`, `rotor`, and `spool`.
# ---------------------------------------------------------------------------


def test_radial_turbine_is_five_streamlines(example_spool):
    spool = example_spool("radial-turbine/radial_turbine-1D.py")["spool"]
    assert spool.num_streamlines == 5


def test_radial_turbine_massflow_attribute_is_the_input_echo(example_spool):
    """As with EEE-HPT: `spool.massflow` is the unconverged 0.1 kg/s guess passed into the
    constructor, not the converged value. See
    `test_eee_hpt_massflow_attribute_is_the_input_echo`.
    """
    spool = example_spool("radial-turbine/radial_turbine-1D.py")["spool"]
    assert spool.massflow == pytest.approx(0.1, rel=RTOL)


# --------------------------------------------------------------------------------------
# RE-BASELINED at PR #27 (`fix/radial-area-and-massflow-sign`).
#
# The three radial-turbine goldens below moved because `compute_streamline_areas` changed:
# the band area is now the exact frustum area pi*(r1 + r2)*slant instead of an expression
# that added a cubic metre to a square metre. See tests/test_current_area_formula.py.
#
# These are CHARACTERIZATION tests. They pin behaviour so that a change cannot pass
# unnoticed; they are not statements that the pinned numbers are correct. They fired exactly
# as designed, and the values are re-taken from the example rather than adjusted to fit.
# Only the radial cases move -- every axial golden in this file is untouched, which is the
# expected signature of a change to the RADIAL branch of the area calculation.
# --------------------------------------------------------------------------------------


def test_radial_turbine_stator_per_streamline(example_spool):
    row = example_spool("radial-turbine/radial_turbine-1D.py")["stator"]
    assert row.P0 == pytest.approx([536657.313099755] * 5, rel=RTOL)
    assert row.T0 == pytest.approx([1222.8779357137] * 5, rel=RTOL)
    # was 0.42080703647747797 before the area fix; P0/T0/Yp are unchanged by it.
    assert row.M == pytest.approx([0.49044906439126945] * 5, rel=RTOL)
    assert row.Yp == pytest.approx([0.0] * 5, abs=1e-12)


def test_radial_turbine_rotor_per_streamline(example_spool):
    """The golden that was flagged as most likely to move once the area fix landed. It did.

    `M` here is the absolute Mach number, and it exceeds 1 at the hub streamline
    (M[0] = 1.2320). That is not a choked passage: choking is set by the meridional Mach
    number, which is far lower at a near-radial station, and a swirling flow can carry an
    absolute Mach number above 1 without it. It is also exactly the kind of station where
    the corrected term in `compute_streamline_areas` is large rather than negligible --
    which is why this row moved and the axial rows did not.
    """
    row = example_spool("radial-turbine/radial_turbine-1D.py")["rotor"]
    # pre-fix P0: [443074.15148622065, 437651.6584138988, 428537.84499390324,
    #              418045.0960861983, 409980.2644715815]
    assert row.P0 == pytest.approx(
        [
            434029.9135235289,
            428744.0618755014,
            419931.65779237944,
            409786.61084366404,
            401992.02228410816,
        ],
        rel=RTOL,
    )
    # pre-fix T0: [1161.8429099481775, 1158.73767317465, 1153.2557718494477,
    #              1146.9118258635594, 1142.398479508869]
    assert row.T0 == pytest.approx(
        [
            1154.970831024012,
            1151.9058186779316,
            1146.5434634683434,
            1140.3418075933444,
            1135.9434687109626,
        ],
        rel=RTOL,
    )
    # pre-fix M: [1.269577634397641, 0.9756593824217092, 0.882204528319701,
    #             0.837409532372009, 0.7773593618539902]
    assert row.M == pytest.approx(
        [
            1.2319911973203372,
            0.9525978451296234,
            0.8626853849786268,
            0.8193645118214383,
            0.7607708835335303,
        ],
        rel=RTOL,
    )
    # Yp is set by the loss model, not the area, and is unchanged by the fix.
    assert row.Yp == pytest.approx([0.1358751871641363] * 5, rel=RTOL)


def test_radial_turbine_convergence_tripwire(example_spool):
    """Coarse tripwire on the pressure-balance solve (see the EEE-HPT equivalent
    above for why this is coarse rather than exact: `minimize_scalar`'s own result
    object is internal to `_balance_pressure` and is not exposed on `spool`).

    The iteration count is unchanged at 15 -- the corrected area moved where the solve
    lands, not how hard it was to get there.
    """
    spool = example_spool("radial-turbine/radial_turbine-1D.py")["spool"]
    assert len(spool.convergence_history) == 15
    last = spool.convergence_history[-1]
    # pre-fix: massflow_std 2.7887726216979658e-05, massflow 0.29864834362636516
    assert last["massflow_std"] == pytest.approx(2.457223793289609e-05, rel=1e-6)
    assert last["massflow"] == pytest.approx(0.3357569515451556, rel=1e-6)


@pytest.mark.slow
def test_the_turbine_examples_are_deterministic(example_spool):
    """If the baseline does not reproduce bit for bit, nothing pinned above means anything.
    Checked here on the radial turbine, the cheaper of the two examples to run twice.
    """
    from tests.conftest import EXAMPLES, _run

    a = _run(EXAMPLES / "radial-turbine" / "radial_turbine-1D.py")[
        "spool"
    ].convergence_history[-1]["massflow_std"]
    b = _run(EXAMPLES / "radial-turbine" / "radial_turbine-1D.py")[
        "spool"
    ].convergence_history[-1]["massflow_std"]
    assert a == b
