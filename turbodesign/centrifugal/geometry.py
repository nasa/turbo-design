"""Meridional flow-path geometry -- where things are.

Two rules, per ``docs/centrifugal/01-design.md`` S3.1 and ``docs/PHYSICS-RULES.md`` rule 4:

1. Parameterise by meridional arc length ``m``, never by ``x``. Differentiating
   ``x(m)`` and ``r(m)`` with respect to ``m`` is well-conditioned everywhere,
   including through a 90-degree bend -- unlike ``dr/dx``, which is singular at a
   radial exit.
2. A :class:`Station` is a quasi-orthogonal cut with an explicit area,
   ``area = 2*pi*r_mean*b*(1-blockage)``. One formula, valid at any meridional
   inclination ``phi``. Never branch on whether a station "looks axial".

This module deliberately does not import anything from ``turbodesign.passage`` or
``turbodesign.flow_math`` -- see ``docs/centrifugal/00-root-cause-analysis.md`` B and F
for the bugs those modules carry (b2 collapsing to zero; a dimensionally inconsistent
band-area formula that goes negative at a radial exit).
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Union

import numpy as np
import numpy.typing as npt

PathLike = Union[str, Path]


@dataclass(frozen=True)
class Station:
    """A quasi-orthogonal cut across the flow passage.

    ``r_hub``/``x_hub`` and ``r_shroud``/``x_shroud`` are the endpoints of the cut.
    ``b`` is the passage width measured *along the cut* (``hypot(dx, dr)``), which is
    the NASA-quoted blade height at a radial station (``docs/PHYSICS-RULES.md`` rule 4) and the
    ordinary annulus height at an axial one -- the same quantity, one formula.
    """

    r_hub: float
    x_hub: float
    r_shroud: float
    x_shroud: float
    r_mean: float
    phi: float  # meridional inclination, radians: 0 = axial, pi/2 = radial
    b: float  # passage width along the cut, m
    blockage: float = 0.0  # fraction of the geometric area blocked (boundary layers, wakes)

    @property
    def area(self) -> float:
        """area = 2*pi*r_mean*b*(1-blockage) -- valid at any phi (docs/PHYSICS-RULES.md rule 4)."""
        return 2.0 * math.pi * self.r_mean * self.b * (1.0 - self.blockage)

    def with_blockage(self, blockage: float) -> "Station":
        """Return a copy of this station with a different blockage fraction applied.

        Geometry (``from_csv`` / ``station_at_radius``) never bakes in a blockage
        assumption -- that is a flow/loss-model choice, applied by the component that
        owns the station (e.g. :class:`turbodesign.centrifugal.components.Impeller`).
        """
        return replace(self, blockage=blockage)


class _Curve:
    """One meridional curve (hub or shroud), parameterised by cumulative arc length.

    ``m[0] = 0``; ``m`` is STRICTLY increasing (every segment ``hypot(dx, dr) > 0``),
    so it is always safe to interpolate -- and to pass to ``np.gradient(..., self.m)``.
    A duplicated coordinate (zero-length segment) would make ``m`` merely
    non-decreasing, and ``np.gradient`` would divide by that zero spacing, producing
    ``phi = nan`` with nothing louder than a ``RuntimeWarning`` -- a landmine that then
    propagates into NaN Vx/Vr downstream. Zero-length segments are dropped (the
    duplicated point is redundant: same (x, r), so no shape information is lost).
    """

    def __init__(self, x: npt.NDArray[np.float64], r: npt.NDArray[np.float64]) -> None:
        x = np.asarray(x, dtype=float)
        r = np.asarray(r, dtype=float)
        if x.shape != r.shape or x.ndim != 1 or x.size < 2:
            raise ValueError("hub/shroud curves need matching 1-D x, r arrays of length >= 2")
        seg = np.hypot(np.diff(x), np.diff(r))
        keep = np.concatenate(([True], seg > 0.0))
        self.x = x[keep]
        self.r = r[keep]
        if self.x.size < 2:
            raise ValueError(
                "curve degenerates to a single point once zero-length segments are dropped"
            )
        seg = np.hypot(np.diff(self.x), np.diff(self.r))
        self.m = np.concatenate(([0.0], np.cumsum(seg)))
        # phi(m) = atan2(dr/dm, dx/dm) -- never dr/dx, which is singular at phi=pi/2.
        dxdm = np.gradient(self.x, self.m)
        drdm = np.gradient(self.r, self.m)
        self._phi = np.arctan2(drdm, dxdm)

    def _phi_at(self, m: float) -> float:
        return float(np.interp(m, self.m, self._phi))

    def _bracket(
        self, driving: npt.NDArray[np.float64], target: float, pick: str
    ) -> tuple[int, float]:
        """Find (i, t) such that target lies between driving[i] and driving[i+1].

        ``pick`` selects which crossing to use when several exist ("first" = nearest
        the start of the curve, "last" = nearest the end) -- e.g. the impeller inlet is
        the first crossing of a given x, the impeller exit is the last crossing of a
        given r.
        """
        diff = driving - target
        sign = np.sign(diff)
        sign[sign == 0.0] = 1.0  # treat an exact hit as a zero-width bracket, not a NaN
        crossings = np.where(np.diff(sign) != 0.0)[0]
        if crossings.size == 0:
            lo, hi = float(driving.min()), float(driving.max())
            raise ValueError(f"target {target} outside curve range [{lo}, {hi}]")
        i = int(crossings[-1] if pick == "last" else crossings[0])
        t = (target - driving[i]) / (driving[i + 1] - driving[i])
        return i, t

    def point_at_r(self, r_target: float, pick: str = "last") -> tuple[float, float, float]:
        """Return (x, r, phi) where this curve crosses r = r_target."""
        i, t = self._bracket(self.r, r_target, pick)
        x = self.x[i] + t * (self.x[i + 1] - self.x[i])
        m = self.m[i] + t * (self.m[i + 1] - self.m[i])
        return float(x), float(r_target), self._phi_at(m)

    def point_at_x(self, x_target: float, pick: str = "first") -> tuple[float, float, float]:
        """Return (x, r, phi) where this curve crosses x = x_target."""
        i, t = self._bracket(self.x, x_target, pick)
        r = self.r[i] + t * (self.r[i + 1] - self.r[i])
        m = self.m[i] + t * (self.m[i + 1] - self.m[i])
        return float(x_target), float(r), self._phi_at(m)


def _read_xr_csv(path: PathLike) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Read a two-column ``x_m,r_m`` CSV (NASA HECC Appendix C convention)."""
    xs: list[float] = []
    rs: list[float] = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if [h.strip().lower() for h in header] != ["x_m", "r_m"]:
            raise ValueError(f"{path}: expected header 'x_m,r_m', got {header}")
        for row in reader:
            if not row:
                continue
            xs.append(float(row[0]))
            rs.append(float(row[1]))
    return np.asarray(xs, dtype=float), np.asarray(rs, dtype=float)


def _station_from_points(
    xh: float, rh: float, phih: float, xs: float, rs: float, phis: float
) -> Station:
    r_mean = 0.5 * (rh + rs)
    b = math.hypot(xh - xs, rh - rs)
    phi = 0.5 * (phih + phis)
    return Station(r_hub=rh, x_hub=xh, r_shroud=rs, x_shroud=xs, r_mean=r_mean, phi=phi, b=b)


class MeridionalPath:
    """The hub and shroud meridional curves of a flow path, and the cuts across them.

    Built from two coordinate tables (NASA Appendix C convention: columns ``x_m,r_m``,
    hub and shroud tabulated independently, generally at different arc-length
    sampling). A :class:`Station` is produced on demand by finding where each curve
    crosses a requested radius (or, for the inducer eye, a requested x) -- never by
    assuming hub and shroud share an index.
    """

    def __init__(self, hub: _Curve, shroud: _Curve) -> None:
        self.hub = hub
        self.shroud = shroud

    @classmethod
    def from_csv(cls, hub_csv_path: PathLike, shroud_csv_path: PathLike) -> "MeridionalPath":
        hub_x, hub_r = _read_xr_csv(hub_csv_path)
        shroud_x, shroud_r = _read_xr_csv(shroud_csv_path)
        return cls(_Curve(hub_x, hub_r), _Curve(shroud_x, shroud_r))

    def station_at_radius(self, r: float) -> Station:
        """The cut where hub and shroud both cross radius ``r``.

        Correct at any inclination: at the impeller exit (radial, phi ~ pi/2) hub and
        shroud share the radius and are separated in x by b2 (docs/PHYSICS-RULES.md rule 4); this
        is found by inverting each curve on r independently, so it needs no "is this
        station radial" branch and no shared-index assumption between the two curves.
        """
        xh, rh, phih = self.hub.point_at_r(r, pick="last")
        xs, rs, phis = self.shroud.point_at_r(r, pick="last")
        return _station_from_points(xh, rh, phih, xs, rs, phis)

    def station_at_x(self, x: float) -> Station:
        """The cut where hub and shroud both cross axial location ``x`` (first crossing).

        The mirror of :meth:`station_at_radius`: at an axial cut hub and shroud share
        ``x`` and differ in ``r`` (``docs/PHYSICS-RULES.md`` rule 4's counterpart at the inlet).
        """
        xh, rh, phih = self.hub.point_at_x(x, pick="first")
        xs, rs, phis = self.shroud.point_at_x(x, pick="first")
        return _station_from_points(xh, rh, phih, xs, rs, phis)

    def inducer_eye(self) -> Station:
        """The impeller inlet station, INFERRED as the axial cut at the hub's minimum radius.

        This is a MODELLING ASSUMPTION, not a geometric fact: it presumes the blade
        leading edge sits exactly at the hub's throat, which is benign in slice 1 (the
        lossless exit state is algebraically independent of the LE station -- PR is
        unchanged for any LE choice) but becomes LOAD-BEARING from slice 4 on, where
        the incidence loss and the inducer-shroud relative Mach number key entirely off
        this station. Callers that know the true LE location should pass
        ``Impeller(x_le=...)`` (see ``turbodesign.centrifugal.components.Impeller``) and
        use :meth:`station_at_x` directly; this heuristic is only the documented
        fallback when no explicit LE is given.

        The hub radius contracts from the inlet duct down to the inducer eye and then
        opens back out through the impeller -- a purely geometric feature, the same on
        any centrifugal flow path, requiring no hand-picked x location. The
        corresponding shroud point is the one sharing that x (the eye is an axial cut:
        hub and shroud share x, differ in r -- the mirror of the radial-exit rule).
        """
        i_eye = int(np.argmin(self.hub.r))
        x_eye = float(self.hub.x[i_eye])
        r_hub = float(self.hub.r[i_eye])
        phi_hub = self.hub._phi_at(float(self.hub.m[i_eye]))
        _, r_shroud, phi_shroud = self.shroud.point_at_x(x_eye, pick="first")
        return _station_from_points(x_eye, r_hub, phi_hub, x_eye, r_shroud, phi_shroud)
