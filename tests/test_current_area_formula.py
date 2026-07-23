# tests/test_current_area_formula.py
"""Tie the reproduction in test_band_area_oracle.py to the library it reproduces.

test_band_area_oracle.py deliberately does not import turbodesign: it is an oracle, and it
has to be able to disagree with the library. The cost of that is a copy of flow_math.py's
band-area expression sitting in the test suite with nothing holding the two together. Once
the formula in the library changes, the copy becomes a stale reproduction of a formula that
no longer exists, and the oracle's tests keep passing while asserting nothing about the code.

This file is the tie. It imports turbodesign and checks that
`flow_math.compute_streamline_areas` still computes what `buggy_band_area` reproduces.

It is expected to fail when the area formula is fixed, and it should then be updated: at
that point `buggy_band_area` is a record of the old formula, `compute_streamline_areas`
should agree with `pappus_band_area` instead, and this file should say so.
"""

import numpy as np
import pytest

from tests.test_band_area_oracle import buggy_band_area, pappus_band_area
from turbodesign.bladerow import BladeRow
from turbodesign.flow_math import compute_streamline_areas

# A cone: one band, hub 0.25 -> shroud 0.35 in radius over 0.10 in x. dx is well above the
# 1e-5 guard at flow_math.py:32, so this takes the radial branch.
X1, R1, X2, R2 = 0.0, 0.25, 0.10, 0.35


def _cone_row() -> BladeRow:
    row = BladeRow()
    row.percent_hub_shroud = np.array([0.0, 1.0])
    row.x = np.array([X1, X2])
    row.r = np.array([R1, R2])
    return row


def test_the_library_still_computes_what_the_oracle_reproduces():
    """If this fails, either the library's area formula changed or the reproduction drifted.
    Either way test_band_area_oracle.py needs updating -- that is what this test is for.
    """
    total_area, streamline_area = compute_streamline_areas(_cone_row())
    expected = buggy_band_area(X1, R1, X2, R2)
    assert streamline_area[1] == pytest.approx(expected, rel=1e-12)
    assert total_area == pytest.approx(expected, rel=1e-12)


def test_the_library_area_is_currently_below_the_geometric_area_on_a_cone():
    """The consequence, stated against the library rather than against the reproduction: on
    this cone the area the library computes is 15.0% below the surface of revolution the
    geometry defines. Fixing flow_math.py:40 will make these two agree and this test fail.
    """
    total_area, _ = compute_streamline_areas(_cone_row())
    exact = pappus_band_area(X1, R1, X2, R2)
    err = (total_area - exact) / exact
    assert err == pytest.approx(-0.150, abs=0.001), (
        f"library area {total_area:.6f} m2 vs geometric area {exact:.6f} m2 ({err:.2%})"
    )
