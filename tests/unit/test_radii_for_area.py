"""Tests for flow_math.radii_for_area - the area-to-passage-geometry bridge used to
size a component's inlet from the non-dimensional mass flow function."""

import numpy as np
import pytest

from turbodesign.flow_math import radii_for_area


@pytest.mark.parametrize("area,mean_radius", [(0.1, 0.3), (0.05, 0.15), (1.2, 0.8)])
def test_recovers_the_requested_area(area, mean_radius):
    r_hub, r_shroud = radii_for_area(area, mean_radius)
    recovered = np.pi * (r_shroud**2 - r_hub**2)
    assert recovered == pytest.approx(area, rel=1e-9)


def test_hub_and_shroud_straddle_the_mean_radius_symmetrically():
    r_hub, r_shroud = radii_for_area(0.08, 0.25)
    assert r_hub < 0.25 < r_shroud
    assert (0.25 - r_hub) == pytest.approx(r_shroud - 0.25, rel=1e-12)
