"""Geometry-module fix pins (physics review of turbodesign/centrifugal, Slice 1).

FIX 5: a duplicated hub/shroud coordinate must not produce a NaN meridional
inclination. FIX 6 (geometry half): ``MeridionalPath.station_at_x`` is the primitive
``Impeller(x_le=...)`` needs, and ``inducer_eye()`` must stay usable as a fallback.
"""

import math

import numpy as np
import pytest

from turbodesign.centrifugal.geometry import MeridionalPath, _Curve


def test_duplicated_coordinate_does_not_produce_nan_phi():
    """FIX 5: a zero-length segment used to divide-by-zero inside np.gradient.

    ``x=[0,1,1,2], r=[1,1,1,1.2]`` has a repeated point at index 1/2 (zero-length
    segment). The old code fed this straight into ``np.gradient(self.x, self.m)``;
    with a zero spacing in ``self.m`` that raises only a ``RuntimeWarning`` and
    returns ``phi = nan`` -- which then poisons Vx/Vr downstream with no loud error.
    Dropping the zero-length segment keeps ``m`` strictly increasing and ``phi``
    finite everywhere.
    """
    curve = _Curve(x=np.array([0.0, 1.0, 1.0, 2.0]), r=np.array([1.0, 1.0, 1.0, 1.2]))
    assert np.all(np.isfinite(curve._phi))
    assert np.all(np.diff(curve.m) > 0.0), "arc length must be STRICTLY increasing"
    # the duplicated point was dropped, not kept as a zero-width step
    assert curve.x.size == 3


def test_curve_rejects_degenerating_to_a_single_point():
    """All points identical: nothing survives dropping zero-length segments."""
    with pytest.raises(ValueError):
        _Curve(x=np.array([1.0, 1.0, 1.0]), r=np.array([2.0, 2.0, 2.0]))


def test_curve_with_no_duplicates_is_unaffected():
    """The common case: no zero-length segments, nothing gets dropped."""
    curve = _Curve(x=np.array([0.0, 1.0, 2.0]), r=np.array([0.5, 0.6, 0.8]))
    assert curve.x.size == 3
    assert np.all(np.isfinite(curve._phi))


# ---------------------------------------------------------------- station_at_x (FIX 6)


@pytest.fixture(scope="module")
def path(data_dir):
    return MeridionalPath.from_csv(
        data_dir / "hecc" / "flowpath_hub.csv",
        data_dir / "hecc" / "flowpath_shroud.csv",
    )


def test_station_at_x_matches_inducer_eye_at_its_own_x(path):
    """station_at_x(x) at the inducer-eye x must reproduce inducer_eye() exactly.

    inducer_eye() is documented as a fallback heuristic (argmin hub radius); an
    explicit x_le is expected to route through station_at_x directly.
    """
    eye = path.inducer_eye()
    explicit = path.station_at_x(eye.x_hub)
    assert explicit.r_hub == pytest.approx(eye.r_hub, rel=1e-9)
    assert explicit.r_shroud == pytest.approx(eye.r_shroud, rel=1e-9)
    assert explicit.area == pytest.approx(eye.area, rel=1e-9)


def test_station_at_x_tracks_a_different_axial_location(path):
    """A station 10 mm upstream of the inducer eye must be a DIFFERENT station."""
    eye = path.inducer_eye()
    upstream = path.station_at_x(eye.x_hub - 0.01)
    assert upstream.x_hub == pytest.approx(eye.x_hub - 0.01, abs=1e-9)
    assert upstream.area != pytest.approx(eye.area, rel=1e-6)


def test_inducer_eye_is_still_the_documented_fallback(path):
    """inducer_eye() keeps working with no x_le supplied -- the fallback is real."""
    eye = path.inducer_eye()
    assert eye.r_hub < eye.r_shroud
    assert math.isfinite(eye.phi)
