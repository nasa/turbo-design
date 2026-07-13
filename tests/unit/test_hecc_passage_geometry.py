"""The HECC meridional passage must not self-intersect.

WHY THIS TEST EXISTS. `nasa-hecc.py` carried fifteen hand-typed endwall
control points that placed the impeller hub and shroud trailing edges at the SAME x
(0.250 m) and separated them in RADIUS by 0.12 mm -- while the comment beside them
claimed ``b2 = 0.01547 m``. That is docs/PHYSICS-RULES.md rule 4 exactly inverted: at a
radial station hub and shroud share the RADIUS and are separated in X.

The consequences were geometric and gross -- an exit area 134x too small, and a passage
whose span went NEGATIVE (hub outboard of shroud) between x = 0.341 and x = 0.361. None
of the 340 tests in this suite could see any of it, because no test looked at the
passage. It was found by a human being LOOKING AT THE PLOT.

That is the visual twin of this project's most-repeated lesson ("a term that cannot fire
cannot be caught being wrong"): a geometry nothing plots is a geometry nothing checks.
This test is the check. It guards the committed flowpath -- the single source of truth
that `turbodesign.centrifugal.geometry.MeridionalPath.from_csv` reads -- against ever
admitting a self-intersecting passage again.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
HECC = REPO / "data" / "hecc"

# NASA/CR-2014-218114/REV1 Table 2. The impeller exit radius, and the passage height
# there -- which is measured ALONG X (rule 4), not along r.
R2_M = 0.21581
B2_M = 0.01547
B2_TOL_M = 0.0005  # 0.5 mm; the published b2 is quoted to ~0.01 mm


def _read(name: str) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    rs: list[float] = []
    with (HECC / name).open() as fh:
        for row in csv.DictReader(fh):
            xs.append(float(row["x_m"]))
            rs.append(float(row["r_m"]))
    return np.asarray(xs), np.asarray(rs)


@pytest.fixture(scope="module")
def flowpath() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    xh, rh = _read("flowpath_hub.csv")
    xs, rs = _read("flowpath_shroud.csv")
    return xh, rh, xs, rs


def test_passage_never_self_intersects(flowpath) -> None:
    """Hub and shroud must stay clear of each other everywhere.

    The defect this catches drove the minimum separation to 0.115 mm and then through
    zero to -0.43 mm. A real passage has a floor of order 10 mm.
    """
    xh, rh, xs, rs = flowpath
    d = np.hypot(xs[None, :] - xh[:, None], rs[None, :] - rh[:, None]).min(axis=1)
    worst = float(d.min())
    where = float(xh[int(np.argmin(d))])
    assert worst > 1e-3, (
        f"passage collapses: minimum hub-to-shroud distance {worst * 1e3:.3f} mm "
        f"at x = {where:.4f} m. A self-intersecting passage has negative area."
    )


def test_b2_is_measured_along_x_not_along_r(flowpath) -> None:
    """docs/PHYSICS-RULES.md rule 4, as a number.

    At the impeller exit both walls sit at r2 and are separated in x by b2. Invert each
    curve on RADIUS independently -- never pair them by index, and never branch on
    "does this station look radial".
    """
    xh, rh, xs, rs = flowpath
    x_hub_at_r2 = float(xh[int(np.argmin(np.abs(rh - R2_M)))])
    x_shroud_at_r2 = float(xs[int(np.argmin(np.abs(rs - R2_M)))])
    b2 = abs(x_hub_at_r2 - x_shroud_at_r2)
    assert b2 == pytest.approx(B2_M, abs=B2_TOL_M), (
        f"b2 measured along x = {b2 * 1e3:.3f} mm, expected {B2_M * 1e3:.3f} mm. "
        "If this reads ~0, hub and shroud have been placed at the same x and "
        "separated in r -- rule 4 inverted, and the exit area will be ~134x too small."
    )


def test_exit_area_matches_the_analytic_value(flowpath) -> None:
    """area = 2*pi*r*b, one formula, valid at any inclination (rule 4)."""
    xh, rh, xs, rs = flowpath
    x_hub_at_r2 = float(xh[int(np.argmin(np.abs(rh - R2_M)))])
    x_shroud_at_r2 = float(xs[int(np.argmin(np.abs(rs - R2_M)))])
    b2 = abs(x_hub_at_r2 - x_shroud_at_r2)
    area = 2.0 * np.pi * R2_M * b2
    assert area == pytest.approx(0.020977, rel=0.01), (
        f"impeller exit area {area:.6f} m^2 against NASA's analytic 0.020977 m^2"
    )
