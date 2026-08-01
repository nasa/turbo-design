# tests/test_current_area_formula.py
"""Tie the reproductions in test_band_area_oracle.py to the library they reproduce.

test_band_area_oracle.py deliberately does not import turbodesign: it is an oracle, and it
has to be able to disagree with the library. The cost of that is a copy of flow_math.py's
band-area expression sitting in the test suite with nothing holding the two together. Once
the formula in the library changes, the copy becomes a stale reproduction of a formula that
no longer exists, and the oracle's tests keep passing while asserting nothing about the code.

This file is the tie.

THE FORMULA HAS NOW BEEN FIXED (PR #27, `fix/radial-area-and-massflow-sign`). This file
previously pinned `compute_streamline_areas` to `buggy_band_area` and stated, in its own
docstring, what to do on the day the fix landed:

    "It is expected to fail when the area formula is fixed, and it should then be updated:
    at that point `buggy_band_area` is a record of the old formula,
    `compute_streamline_areas` should agree with `pappus_band_area` instead, and this file
    should say so."

This is that update, and this file now says so. `compute_streamline_areas` is checked
against `pappus_band_area` -- the exact lateral area of the frustum the two points define,
pi*(r1 + r2)*slant. `buggy_band_area` is retained in the oracle purely as a record of the
old expression, and is used here in the opposite direction: to assert the library has NOT
regressed back to it.

Why the old one was wrong, for the record. It read

    S = r2 - r1                       # a LENGTH, in metres -- not a slope
    C = sqrt(1 + (S/dx)**2)
    area = 2*pi*C*(S/2 * dx**2 + r1*dx)

`S/2 * dx**2` is a cubic metre and it was added to `r1*dx`, a square metre. A sum whose
terms carry different dimensions cannot be right in any unit system, and the consequence is
not small: the expression is exact only for a cylinder (S = 0, where the bad term vanishes),
and on the cone below it was 15% low.
"""

import numpy as np
import pytest

from tests.test_band_area_oracle import (
    buggy_band_area,
    pappus_band_area,
    slope_variant_band_area,
)
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


def test_the_library_computes_the_exact_geometric_area():
    """`compute_streamline_areas` must equal the surface of revolution the geometry defines.

    If this fails, either the library's area formula changed again or the Pappus
    reproduction drifted. Either way test_band_area_oracle.py needs updating -- that is what
    this test is for.
    """
    total_area, streamline_area = compute_streamline_areas(_cone_row())
    exact = pappus_band_area(X1, R1, X2, R2)
    assert streamline_area[1] == pytest.approx(exact, rel=1e-12)
    assert total_area == pytest.approx(exact, rel=1e-12)


def test_the_library_has_not_regressed_to_the_dimensionally_broken_formula():
    """The tripwire in the other direction.

    `buggy_band_area` is 15% low on this cone. Pinning the gap keeps the old expression
    from creeping back in unnoticed -- and keeps `buggy_band_area` itself honest, since a
    reproduction nothing exercises is a reproduction nobody notices going stale.
    """
    total_area, _ = compute_streamline_areas(_cone_row())
    old = buggy_band_area(X1, R1, X2, R2)
    err = (old - total_area) / total_area
    assert err == pytest.approx(-0.150, abs=0.001), (
        f"the superseded formula gives {old:.6f} m2 against the library's {total_area:.6f} m2 "
        f"({err:.2%}). If this gap has closed, the library may have regressed to it."
    )


def test_the_area_is_dimensionally_homogeneous():
    """Scale the whole geometry by k; a true area must scale by k**2, exactly.

    This is the check that would have caught the old formula on day one, without knowing
    anything about frustums: it summed a cubic metre and a square metre, so it scaled as
    neither. It came out ~5.97e6 across a 1e3 length scaling where an area must give 1e6.
    """
    k = 1000.0
    plain, _ = compute_streamline_areas(_cone_row())

    scaled_row = BladeRow()
    scaled_row.percent_hub_shroud = np.array([0.0, 1.0])
    scaled_row.x = np.array([X1 * k, X2 * k])
    scaled_row.r = np.array([R1 * k, R2 * k])
    scaled, _ = compute_streamline_areas(scaled_row)

    assert scaled / plain == pytest.approx(k**2, rel=1e-12)


def test_the_library_returns_a_SIGNED_area_not_a_geometric_one():
    """The library's area is signed: reverse the cut and it changes sign.

    This is deliberate in PR #27 -- `# Signed area: sign follows dx to maintain massflow
    sign convention` -- and this test does not argue with it. It pins it, because nothing
    else does, and because the sign is invisible in every forward-ordered case the rest of
    this file exercises.

    Note what the suite already says about this shape. `slope_variant_band_area` in the
    oracle is "the obvious repair to buggy_band_area ... It fixes the dimensional error but
    not the sign", and `test_area_is_never_negative_for_a_reversed_cut` records the verdict:
    a reversed cut "returns a negative area, which Pappus -- and an area -- cannot."

    The library now computes that variant EXACTLY, reversed cuts included. So the geometric
    magnitude is right and agrees with Pappus, while the sign is a flow-direction convention
    carried on the same number as the area. That is a real design decision and it may well be
    the intended one; it is pinned here so that changing it cannot pass unnoticed, and so the
    disagreement with the geometric oracle is stated rather than latent.
    """
    forward, _ = compute_streamline_areas(_cone_row())

    reversed_row = BladeRow()
    reversed_row.percent_hub_shroud = np.array([0.0, 1.0])
    reversed_row.x = np.array([X2, X1])
    reversed_row.r = np.array([R2, R1])
    backward, _ = compute_streamline_areas(reversed_row)

    # Magnitude: agrees with the geometric oracle in both directions.
    exact = pappus_band_area(X1, R1, X2, R2)
    assert abs(forward) == pytest.approx(exact, rel=1e-12)
    assert abs(backward) == pytest.approx(exact, rel=1e-12)

    # Sign: the convention. Pappus is unsigned and cannot express this.
    assert forward > 0.0
    assert backward < 0.0
    assert backward == pytest.approx(-forward, rel=1e-12)

    # And it is precisely the variant the oracle names and rejects as a geometric area.
    assert forward == pytest.approx(slope_variant_band_area(X1, R1, X2, R2), rel=1e-12)
    assert backward == pytest.approx(slope_variant_band_area(X2, R2, X1, R1), rel=1e-12)
