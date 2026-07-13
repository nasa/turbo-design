# tests/test_axial_characterization.py
"""Pin the current behaviour of the axial compressor examples.

These are characterization goldens: they record what the code does today, not what is
physically correct. Several are expected to move when known defects are fixed; the point of
pinning them is that the change is then visible and has to be stated.

Pinned per streamline (hub, mean, shroud). Both examples run num_streamlines = 3, the
CompressorSpool default. A meanline-only golden would be blind to any defect whose error
vanishes at constant radius -- and both examples hold the mean radius constant.
"""

import numpy as np
import pytest

RTOL = 1e-9


def test_example_9_1_is_three_streamlines(example_spool):
    spool = example_spool("mattingly-axial-compressor/example9.1.py")["spool"]
    assert spool.num_streamlines == 3


def test_example_9_1_overall(example_spool):
    spool = example_spool("mattingly-axial-compressor/example9.1.py")["spool"]
    assert spool.massflow == pytest.approx(22.6796, rel=RTOL)
    assert spool.overall_pressure_ratio() == pytest.approx(1.3360587343032544, rel=RTOL)


def test_example_9_1_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.1.py")["rotor"]
    assert rotor.r == pytest.approx([0.274612394592, 0.3048, 0.334987605408], rel=RTOL)
    assert rotor.P0 == pytest.approx(
        [118723.70502439793, 135158.28238870634, 147801.11230809323], rel=RTOL
    )
    assert rotor.T0 == pytest.approx(
        [302.717182747679, 312.798819587163, 319.179654224343], rel=RTOL
    )
    assert rotor.P0R == pytest.approx(
        [100109.87642584258, 106221.6927355795, 105585.90931530058], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.256792735336, 291.899539873286, 289.806659750478], rel=RTOL
    )
    assert rotor.M == pytest.approx(
        [0.82329899123, 0.882171122458, 0.947232789356], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.638851397703, 0.627935526295, 0.597910853155], rel=RTOL
    )


def test_example_9_1_stator_per_streamline(example_spool):
    stator = example_spool("mattingly-axial-compressor/example9.1.py")["stator"]
    assert stator.P0 == pytest.approx(
        [122755.33113047821, 135108.6406580774, 148376.59867602392], rel=RTOL
    )
    assert stator.M == pytest.approx(
        [0.586561267717, 0.670423034066, 0.740061151736], rel=RTOL
    )


def test_example_9_2_overall(example_spool):
    spool = example_spool("mattingly-axial-compressor/example9.2.py")["spool"]
    assert spool.massflow == pytest.approx(22.68, rel=RTOL)
    assert spool.overall_pressure_ratio() == pytest.approx(1.292878358443007, rel=RTOL)


def test_example_9_2_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.2.py")["rotor"]
    assert rotor.r == pytest.approx([0.273978687961, 0.3048, 0.335621312039], rel=RTOL)
    assert rotor.P0R == pytest.approx(
        [97543.68041388667, 103111.38595056451, 102563.81089289424], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.045210087425, 291.849742618825, 289.967113368956], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.650202734506, 0.636853067701, 0.605750066296], rel=RTOL
    )


def test_the_annulus_contracts_off_the_meanline(example_spool):
    """What the per-streamline pin above is for.

    Both examples hold the mean radius constant while contracting the annulus, so r2/r1 is
    exactly 1.000 at the mean and differs from 1 at hub and shroud. Any error that scales
    with (U2 - U1) is therefore absent from the meanline and present only off it.
    """
    for name in (
        "mattingly-axial-compressor/example9.1.py",
        "mattingly-axial-compressor/example9.2.py",
    ):
        loc = example_spool(name)
        inlet, rotor = loc["spool"].blade_rows[0], loc["rotor"]
        ratio = np.asarray(rotor.r, float) / np.asarray(inlet.r, float)
        assert ratio[1] == pytest.approx(1.0, abs=1e-12), f"{name}: mean radius moved"
        assert abs(ratio[0] - 1.0) > 1e-3, f"{name}: hub does not contract"
        assert abs(ratio[2] - 1.0) > 1e-3, f"{name}: shroud does not contract"


def test_the_examples_are_deterministic(example_spool):
    """If the baseline does not reproduce bit for bit, nothing pinned above means anything."""
    from tests.conftest import EXAMPLES, _run

    a = _run(EXAMPLES / "mattingly-axial-compressor" / "example9.1.py")[
        "spool"
    ].overall_pressure_ratio()
    b = _run(EXAMPLES / "mattingly-axial-compressor" / "example9.1.py")[
        "spool"
    ].overall_pressure_ratio()
    assert a == b
