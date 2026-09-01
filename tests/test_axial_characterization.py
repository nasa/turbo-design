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
    # Re-baselined: stator_calc previously never re-derived a stator's P0 from its
    # upstream rotor's converged P0 (only from initialize()'s one-time, pre-convergence
    # seed), so overall_pressure_ratio() - which reads the last row's P0 - was wrong
    # whenever that last row is a stator, as it is here.
    spool = example_spool("mattingly-axial-compressor/example9.1.py")["spool"]
    assert spool.massflow == pytest.approx(22.6796, rel=RTOL)
    assert spool.overall_pressure_ratio() == pytest.approx(1.322718069014324, rel=RTOL)


def test_example_9_1_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.1.py")["rotor"]
    assert rotor.r == pytest.approx([0.274612394592, 0.3048, 0.334987605408], rel=RTOL)
    assert rotor.P0 == pytest.approx(
        [119628.40319595579, 134654.94917312535, 145123.17232950567], rel=RTOL
    )
    assert rotor.T0 == pytest.approx(
        [302.20998821723526, 312.6535136186777, 319.44470528076056], rel=RTOL
    )
    assert rotor.P0R == pytest.approx(
        [101472.66241132512, 106007.04371937654, 103382.95293928268], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.26371192569735, 291.9065465214719, 289.81361672926107], rel=RTOL
    )
    assert rotor.M == pytest.approx(
        [0.8300590707928704, 0.8836658122756325, 0.9454460873879358], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.654180180595538, 0.6319804100205337, 0.5917381359155736], rel=RTOL
    )


def test_example_9_1_stator_per_streamline(example_spool):
    # Re-baselined along with test_example_9_1_overall - see that test's comment.
    stator = example_spool("mattingly-axial-compressor/example9.1.py")["stator"]
    assert stator.P0 == pytest.approx(
        [121107.41535333314, 134654.94917312535, 146421.85771165002], rel=RTOL
    )
    assert stator.M == pytest.approx(
        [0.5963208015278048, 0.6888803194859718, 0.7469872531703152], rel=RTOL
    )


def test_example_9_2_overall(example_spool):
    # Re-baselined: see test_example_9_1_overall's comment - same stator P0 fix.
    spool = example_spool("mattingly-axial-compressor/example9.2.py")["spool"]
    assert spool.massflow == pytest.approx(22.68, rel=RTOL)
    assert spool.overall_pressure_ratio() == pytest.approx(1.265613361918456, rel=RTOL)


def test_example_9_2_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.2.py")["rotor"]
    assert rotor.r == pytest.approx([0.273978687961, 0.3048, 0.335621312039], rel=RTOL)
    assert rotor.P0R == pytest.approx(
        [98740.7993354694, 102918.7847158074, 100612.12365309549], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.0527983355748, 291.85741965922836, 289.97441139016513], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.6644070309683578, 0.6407609185183611, 0.6002948269776887], rel=RTOL
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
