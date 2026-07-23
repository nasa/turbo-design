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
    assert spool.overall_pressure_ratio() == pytest.approx(1.3341973512851182, rel=RTOL)


def test_example_9_1_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.1.py")["rotor"]
    assert rotor.r == pytest.approx([0.274612394592, 0.3048, 0.334987605408], rel=RTOL)
    assert rotor.P0 == pytest.approx(
        [119628.39581103637, 134654.9438487545, 145123.1800222211], rel=RTOL
    )
    assert rotor.T0 == pytest.approx(
        [302.2099828302515, 312.6535100310146, 319.4447100979148], rel=RTOL
    )
    assert rotor.P0R == pytest.approx(
        [101472.6344434061, 106007.01485885675, 103382.93537587197], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.26368910588, 291.906523697345, 289.8136025909995], rel=RTOL
    )
    assert rotor.M == pytest.approx(
        [0.8300588723044268, 0.8836656583201487, 0.9454459831924451], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.6541796871156736, 0.6319799253120101, 0.5917376972636], rel=RTOL
    )


def test_example_9_1_stator_per_streamline(example_spool):
    stator = example_spool("mattingly-axial-compressor/example9.1.py")["stator"]
    assert stator.P0 == pytest.approx(
        [124686.40724035529, 135017.66232236876, 145970.53079910018], rel=RTOL
    )
    assert stator.M == pytest.approx(
        [0.6092420162196824, 0.6720905746225406, 0.7269550972159657], rel=RTOL
    )


def test_example_9_2_overall(example_spool):
    spool = example_spool("mattingly-axial-compressor/example9.2.py")["spool"]
    assert spool.massflow == pytest.approx(22.68, rel=RTOL)
    assert spool.overall_pressure_ratio() == pytest.approx(1.2911654599889533, rel=RTOL)


def test_example_9_2_rotor_per_streamline(example_spool):
    rotor = example_spool("mattingly-axial-compressor/example9.2.py")["rotor"]
    assert rotor.r == pytest.approx([0.273978687961, 0.3048, 0.335621312039], rel=RTOL)
    assert rotor.P0R == pytest.approx(
        [98740.66209918521, 102918.64302259442, 100612.03719717609], rel=RTOL
    )
    assert rotor.T0R == pytest.approx(
        [288.05267532349905, 291.8572966174646, 289.97433510707276], rel=RTOL
    )
    assert rotor.M_rel == pytest.approx(
        [0.6644041789821085, 0.6407581013471035, 0.6002922414935745], rel=RTOL
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
