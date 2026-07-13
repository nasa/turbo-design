"""Tests for my_scripts/flowpath_builder.py (WAVE B geometry reconstruction).

Covers: the superellipse turn's endpoint tangents (exact, analytic -- see the module
docstring's derivation), the b2 sign convention (docs/PHYSICS-RULES.md rule 4), a CSV
round trip through the real geometry entry point (turbodesign.centrifugal.geometry --
NOT modified by this test, only exercised), the constant-b vaneless diffuser, and the
analytic quarter-circle check that anchors ``meridional_length``'s discretisation error.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "my_scripts"))

from flowpath_builder import (  # noqa: E402
    _superellipse_turn,
    _tangent_slopes,
    build_flowpath,
    eckardt_rotor_o_lm_check,
    meridional_length,
    write_flowpath_csv,
)

from turbodesign.centrifugal.geometry import MeridionalPath  # noqa: E402

# A representative meanline table, loosely HECC-scaled, used across several tests.
_R1H, _R1S = 0.04, 0.11
_R2, _B2 = 0.15, 0.012
_LZ = 0.06
_R_EXIT = 0.20
_INLET_DUCT_LEN = 0.05


# --------------------------------------------------------------------------- tangents


@pytest.mark.parametrize("shape", [1.5, 2.0, 3.0])
def test_le_tangent_is_axial_and_te_tangent_is_radial(shape):
    """dr/dx == 0.0 exactly at the LE; dx/dr == 0.0 exactly at the TE -- for any shape > 1.

    These are exact IEEE values from the closed-form derivative (module docstring), not
    numerical approximations, so they are asserted with ``==`` rather than ``approx``.
    """
    a, b = 0.08, 0.05
    dr_dx_le, dx_dr_le = _tangent_slopes(0.0, a, b, shape)
    assert dr_dx_le == 0.0
    assert math.isinf(dx_dr_le)

    dr_dx_te, dx_dr_te = _tangent_slopes(1.0, a, b, shape)
    assert math.isinf(dr_dx_te)
    assert dx_dr_te == 0.0


@pytest.mark.parametrize("shape", [1.5, 2.0, 3.0])
def test_sampled_curve_endpoints_match_le_te(shape):
    """The sampled curve's first/last points reproduce (x_le, r_le) / (x_te, r_te).

    Default (relative) ``pytest.approx`` tolerance, not ``abs=0.0``: r_le=0.02 is computed
    as ``r_te - B`` with ``B = r_te - r_le`` -- an algebraic identity, but not bit-exact
    under IEEE float subtraction, so a handful of ULPs of difference is expected and fine.
    """
    x, r = _superellipse_turn(0.0, 0.02, 0.10, 0.20, shape, 50)
    assert x[0] == pytest.approx(0.0, abs=1e-15)
    assert r[0] == pytest.approx(0.02)
    assert x[-1] == pytest.approx(0.10)
    assert r[-1] == pytest.approx(0.20)


# --------------------------------------------------------------------------- b2 / rule 4


def test_hub_and_shroud_share_radius_and_are_separated_by_exactly_b2_at_te():
    """docs/PHYSICS-RULES.md rule 4: hub/shroud share r2, separated ALONG X by b2.

    Sign convention (module docstring): hub sits DOWNSTREAM of shroud at the TE, i.e.
    x_hub_TE - x_shroud_TE == +b2. Confirmed against data/hecc/flowpath_{hub,shroud}.csv.
    """
    hub_xr, shroud_xr = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, r_exit=_R2, inlet_duct_len=0.0, lz_refers_to="hub"
    )
    # last sample of each curve is the TE (r_exit == r2 here, so no diffuser was appended).
    x_hub_te, r_hub_te = hub_xr[-1]
    x_shroud_te, r_shroud_te = shroud_xr[-1]
    assert r_hub_te == pytest.approx(_R2, abs=1e-12)
    assert r_shroud_te == pytest.approx(_R2, abs=1e-12)
    assert r_hub_te == pytest.approx(r_shroud_te, abs=1e-12)
    assert (x_hub_te - x_shroud_te) == pytest.approx(_B2, abs=1e-12)


# --------------------------------------------------------------------------- CSV round trip


def test_csv_round_trip_reproduces_b2_and_area_at_r2(tmp_path):
    """The emitted CSV, read back through the REAL geometry entry point
    (MeridionalPath.from_csv), gives station_at_radius(r2).b == b2 and
    .area == 2*pi*r2*b2 -- docs/PHYSICS-RULES.md rule 4: one area formula, valid at any phi.
    """
    hub_xr, shroud_xr = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, _R_EXIT, _INLET_DUCT_LEN, lz_refers_to="hub"
    )
    hub_path, shroud_path = write_flowpath_csv("case", hub_xr, shroud_xr, out_dir=tmp_path)

    path = MeridionalPath.from_csv(hub_path, shroud_path)
    station = path.station_at_radius(_R2)

    assert station.b == pytest.approx(_B2, abs=1e-9)
    assert station.area == pytest.approx(2.0 * math.pi * _R2 * _B2, rel=1e-9)


# --------------------------------------------------------------------------- monotonicity


def test_r_increases_monotonically_along_both_curves():
    """A meanline solver walks the flow path by radius; r must never go backwards.

    Non-decreasing, not strictly increasing: the inlet-duct point and the impeller curve's
    own t=0 sample share the same r (both are the LE radius -- the duct is flat and the
    impeller's tangent is exactly axial there too, docstring derivation), so that one join
    has dr=0. Everywhere else r is strictly increasing.
    """
    hub_xr, shroud_xr = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, _R_EXIT, _INLET_DUCT_LEN, lz_refers_to="hub"
    )
    assert np.all(np.diff(hub_xr[:, 1]) >= 0.0)
    assert np.all(np.diff(shroud_xr[:, 1]) >= 0.0)
    assert hub_xr[-1, 1] > hub_xr[0, 1]
    assert shroud_xr[-1, 1] > shroud_xr[0, 1]


# --------------------------------------------------------------------------- diffuser


def test_diffuser_walls_are_parallel_constant_b(tmp_path):
    """b(r) is exactly b2 from r2 out to r_exit -- the vaneless diffuser has parallel walls."""
    hub_xr, shroud_xr = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, _R_EXIT, inlet_duct_len=0.0, lz_refers_to="hub"
    )
    hub_path, shroud_path = write_flowpath_csv("diffuser-case", hub_xr, shroud_xr, out_dir=tmp_path)
    path = MeridionalPath.from_csv(hub_path, shroud_path)
    for r in np.linspace(_R2, _R_EXIT, 6):
        station = path.station_at_radius(float(r))
        assert station.b == pytest.approx(_B2, abs=1e-9), f"b(r={r}) != b2"


# --------------------------------------------------------------------------- Lz ambiguity


def test_lz_refers_to_hub_vs_shroud_changes_te_x_by_b2():
    """Switching lz_refers_to shifts which curve's TE plane equals Lz -- by exactly b2."""
    hub_h, shroud_h = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, r_exit=_R2, inlet_duct_len=0.0, lz_refers_to="hub"
    )
    hub_s, shroud_s = build_flowpath(
        _R1H, _R1S, _R2, _B2, _LZ, r_exit=_R2, inlet_duct_len=0.0, lz_refers_to="shroud"
    )
    assert hub_h[-1, 0] == pytest.approx(_LZ, abs=1e-12)
    assert shroud_s[-1, 0] == pytest.approx(_LZ, abs=1e-12)
    assert hub_h[-1, 0] != pytest.approx(hub_s[-1, 0], abs=1e-9)
    assert shroud_h[-1, 0] != pytest.approx(shroud_s[-1, 0], abs=1e-9)


# --------------------------------------------------------------------------- meridional_length


def test_meridional_length_reproduces_quarter_circle_for_shape_2_A_equals_B():
    """shape=2, A=B is a TRUE circular quadrant (algebraically: x^2+(R-r)^2=R^2), so the
    mean-line arc length must reproduce the analytic quarter-circle length pi/2 * R.

    Tolerance: meridional_length integrates piecewise-linear CHORDS between sampled points,
    which underestimates a convex arc's true length; the error scales like O(1/n_points^2)
    (a standard trapezoidal/chord discretisation result for a smooth curve). Measured
    empirically at n_points=200 (this test's setting) the relative error is ~4e-5; rel=1e-3
    is used here as a comfortable (25x) margin above that measurement, not a tuned value.
    """
    radius = 0.1
    n_points = 200
    curve_x, curve_r = _superellipse_turn(0.0, 0.0, radius, radius, 2.0, n_points)
    curve_xr = np.column_stack([curve_x, curve_r])

    length = meridional_length(curve_xr, curve_xr)
    exact = 0.5 * math.pi * radius
    assert length == pytest.approx(exact, rel=1e-3)


def test_meridional_length_mean_line_of_two_different_curves_is_between_them():
    """Sanity check on the pairing-by-arc-fraction: the mean line's length sits between the
    (very different) hub and shroud impeller-turn lengths for an asymmetric case."""
    hub_x, hub_r = _superellipse_turn(0.0, 0.0, 0.20, 0.20, 2.0, 100)
    shroud_x, shroud_r = _superellipse_turn(0.0, 0.0, 0.05, 0.05, 2.0, 100)
    hub_xr = np.column_stack([hub_x, hub_r])
    shroud_xr = np.column_stack([shroud_x, shroud_r])

    hub_len = meridional_length(hub_xr, hub_xr)
    shroud_len = meridional_length(shroud_xr, shroud_xr)
    mean_len = meridional_length(hub_xr, shroud_xr)

    assert shroud_len < mean_len < hub_len


# --------------------------------------------------------------------------- Eckardt rotor O check


def test_eckardt_rotor_o_default_reconstruction_is_within_a_few_percent_of_published_lm():
    """Kovář et al. 2021 Table 1 publishes Lm = 0.1739 m for Eckardt rotor O. The DEFAULT
    reconstruction (shape=2, lz_refers_to='hub') must land within a few percent of that --
    this is the independent check the module exists to provide, not a tuned match."""
    report = eckardt_rotor_o_lm_check()
    assert report["Lm_published_m"] == pytest.approx(0.1739)
    assert abs(report["disagreement_pct"]) < 2.0
    # the shape that WOULD match exactly is reported too, and should be close to (not
    # forced to equal) the untouched default of 2.0 -- evidence the default is reasonable,
    # not a demonstration that it is exact.
    assert 1.5 < report["shape_matching_published_Lm"] < 2.5


# --------------------------------------------------------------------------- input validation


def test_shape_must_exceed_one():
    with pytest.raises(ValueError):
        build_flowpath(_R1H, _R1S, _R2, _B2, _LZ, _R2, 0.0, shape=1.0)


def test_r_exit_below_r2_is_rejected():
    with pytest.raises(ValueError):
        build_flowpath(_R1H, _R1S, _R2, _B2, _LZ, r_exit=_R2 - 0.01, inlet_duct_len=0.0)


def test_invalid_lz_refers_to_is_rejected():
    with pytest.raises(ValueError):
        build_flowpath(_R1H, _R1S, _R2, _B2, _LZ, _R2, 0.0, lz_refers_to="nose")
