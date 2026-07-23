"""Reconstruct a hub/shroud meridional flow-path CSV pair from a published meanline table.

WHY THIS EXISTS
----------------
``turbodesign.centrifugal.geometry.MeridionalPath.from_csv(hub_csv, shroud_csv)`` is the
ONLY geometry entry point in ``turbodesign/centrifugal/`` (see that module's docstring and
``docs/PHYSICS-RULES.md`` rule 4) and it stays that way -- this script does not import from,
or modify, anything under ``turbodesign/``. ``from_csv`` reads two-column ``x_m,r_m`` CSVs
(SI metres, NASA HECC Appendix C convention).

Most published centrifugal-impeller meanline tables give a handful of scalars -- inducer
hub/shroud radius, exit radius, exit passage width, axial length -- not the coordinate
tables ``from_csv`` needs (Kovář et al. 2021, *Energies* 14(24) 8545, Table 1, is exactly
this: it tabulates ``D1h, D1s, D2, b2, Lz, Lm`` for the Eckardt rotors, no coordinates).
This module bridges that gap: it RECONSTRUCTS a plausible hub and shroud curve from the
scalar table, using a shape parameter that is explicitly a free knob, not a hidden
assumption -- and it ships an independent check of how good that reconstruction is,
against a number the source paper itself publishes (see ``eckardt_rotor_o_lm_check``
below). Agreeing with the paper's own published meridional length is evidence; agreeing
with our own estimate at a hand-picked shape is not -- that distinction is the entire
methodology of this project (docs/PHYSICS-RULES.md, "Validation doctrine").

SIGN CONVENTION -- confirmed by reading data/hecc/flowpath_hub.csv and flowpath_shroud.csv
--------------------------------------------------------------------------------------------
docs/PHYSICS-RULES.md rule 4: at a RADIAL station (impeller exit) hub and shroud share the
same radius and are separated in x by b2. Which one sits downstream is not stated by the
rule -- it has to come from data. It was determined here by finding, in each CSV, the
run of samples where x is (numerically) constant while r sweeps a wide range -- i.e. the
parallel-wall vaneless-diffuser segment immediately downstream of the impeller TE, where
each wall sits at a fixed x and only r changes:

    hub:    x is EXACTLY 0.13376277539999998 m for 5 consecutive rows (176-180), r sweeping
            0.2215 -> 0.2539 m.
    shroud: x is CONSTANT to within +/-2e-4 m (digitisation noise) around 0.120004 m for
            ~30 rows (173-205), r sweeping 0.216 -> 0.278 m.

    x_hub_TE - x_shroud_TE = 0.13376 - 0.12000 = +0.01376 m  ( > 0 )

**The HUB trailing edge sits at LARGER x than the SHROUD trailing edge** -- the hub runs
axially further downstream before completing its turn to radial. This module follows that
sign throughout: ``x_te_hub = x_te_shroud + b2``, always, regardless of which curve the
``Lz`` argument is anchored to (see ``lz_refers_to`` below).

THE SUPERELLIPSE TURN -- derivation
------------------------------------
Hub and shroud each turn from AXIAL at the LE (tangent dr/dx = 0: the meridional direction
is parallel to the x-axis) to RADIAL at the TE (tangent dx/dr = 0: the meridional direction
is parallel to the r-axis). A task brief for this module suggested two candidate closed
forms; neither was accepted on faith -- both are wrong or under-specified, and the one
used here was re-derived and checked analytically (see ``_tangent_slopes`` and
``tests/unit/test_flowpath_builder.py``).

Start from the superellipse quadrant through (x_LE, r_LE) and (x_TE, r_TE), with
``A = x_TE - x_LE`` and ``B = r_TE - r_LE``:

    ((x - x_LE) / A)^s + ((r_TE - r) / B)^s = 1                       (s = ``shape``, s > 1)

Parameterise by ``t`` in [0, 1] via ``(x - x_LE)/A = t`` (i.e. sample x LINEARLY -- this
keeps the parametrisation explicit and lets the CSV point density be controlled directly
by ``n_points``). Substituting and solving for r:

    x(t) = x_LE + A * t
    r(t) = r_TE - B * (1 - t**s) ** (1/s)

Check: t=0 -> x=x_LE, r = r_TE - B*(1)**(1/s) = r_TE - B = r_LE.  (LE, correct)
       t=1 -> x=x_TE, r = r_TE - B*(0)**(1/s) = r_TE.             (TE, correct)

Endpoint tangents, differentiating r(t) (dx/dt = A, constant):

    dr/dt = B * t**(s-1) * (1 - t**s)**(1/s - 1)

    t -> 0:  t**(s-1) -> 0  for s > 1   =>  dr/dt = 0   =>  dr/dx = 0        (axial tangent)
    t -> 1:  (1-t**s)**(1/s-1), with 1/s-1 < 0 for s > 1, base -> 0+
                                        =>  dr/dt -> +inf  =>  dx/dr = 0     (radial tangent)

Both limits are exact IEEE values under floating point (``0.0**positive = 0.0``,
``finite/inf = 0.0``), not numerical approximations -- see ``_tangent_slopes``, asserted
exactly (not `pytest.approx`) in the test suite. This holds for ANY ``shape > 1``; ``shape``
is a free parameter with no "correct" value -- see the module-level discussion of degrees
of freedom below. ``shape = 2`` is the plain ellipse quadrant. Larger ``shape`` pushes the
turn later (a "squarer" corner, approaching the L-shaped bounding-box path of length
``A + B`` as shape -> inf); smaller ``shape`` (-> 1) pulls the turn earlier, approaching the
straight chord of length ``hypot(A, B)`` (the shortest possible connector, and the one
``shape = 1`` degenerates to -- excluded here because it has no flat LE/TE tangent).

DEGREES OF FREEDOM IN THE RECONSTRUCTION -- made explicit, not buried
-----------------------------------------------------------------------
1. ``shape`` -- see above. THE free parameter; it is meant to be perturbed and its
   sensitivity reported, not fixed once and forgotten.
2. ``lz_refers_to`` -- published meanline tables give ONE axial-length number ``Lz`` for
   the whole impeller, but hub and shroud have DIFFERENT axial extents (their TE planes
   are offset by ``b2``, per the sign convention above) and a table essentially never says
   which wall ``Lz`` was measured along. Kovář et al. 2021 Table 1 is exactly this case:
   one ``Lz = 130 mm`` for both Eckardt rotors, no statement of which contour it follows.
   This module makes the choice an explicit argument, ``lz_refers_to in {"hub", "shroud"}``,
   rather than silently picking one.

   **Default: ``"hub"``.** This was NOT read directly off the HECC CSVs -- HECC's own
   ``data/hecc/design_point.csv`` carries no published ``Lz`` or ``Lm`` scalar to check a
   convention against, only the coordinate tables (which give a sign convention -- see
   above -- but not a scalar-table semantic). The default is instead justified by the one
   place in this project's data where BOTH conventions can be checked against an
   independently PUBLISHED meridional length: Kovář et al. 2021 Table 1 for Eckardt rotor O
   publishes ``Lz = 0.130 m`` and, separately, ``Lm = 0.1739 m``. At the plain-ellipse
   default (``shape = 2``, chosen with no knowledge of this check), interpreting
   ``Lz`` as the HUB's axial extent reproduces the published ``Lm`` to **+0.56%**;
   interpreting it as the SHROUD's axial extent misses by **+12.7%**. See
   ``eckardt_rotor_o_lm_check`` for the numbers and the two-sided comparison that produced
   this default -- it was not tuned to match; both conventions were tried at the same
   untouched ``shape = 2`` and one was simply far closer.

GEOMETRY LAYOUT
-----------------
    x_LE = 0 (arbitrary reference; hub and shroud share it -- the mirror of the b2 rule,
              rule 4's "at an AXIAL station hub and shroud share x and differ in r").

    [inlet duct, len=inlet_duct_len]  straight axial line, r = r1h / r1s constant.
    [impeller turn, len from Lz/lz_refers_to]  superellipse quadrant, r1h/r1s -> r2.
    [vaneless diffuser, out to r_exit]  straight RADIAL line at fixed x = x_te_hub / x_te_shroud
              (parallel walls, b = b2 constant -- docs/PHYSICS-RULES.md rule 4's area formula
              is one formula valid at any phi, so a constant-b vaneless diffuser is just this
              same curve-pair machinery with A = 0 for that segment).

Straight segments use 2 points each: a straight line interpolates EXACTLY between its two
endpoints (the geometry module's ``Station`` lookups are linear-interpolation-based -- see
``turbodesign/centrifugal/geometry.py``), so no density is lost by keeping them minimal.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Literal, Tuple

import numpy as np
import numpy.typing as npt

# A curve as returned internally: (x, r), each a 1-D float array of matching length.
_XR = Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]

LzRef = Literal["hub", "shroud"]


# --------------------------------------------------------------------------- superellipse


def _superellipse_turn(
    x_le: float, r_le: float, x_te: float, r_te: float, shape: float, n_points: int
) -> _XR:
    """Sample the axial-to-radial superellipse turn described in the module docstring.

    ``x(t) = x_le + A*t``, ``r(t) = r_te - B*(1 - t**shape)**(1/shape)``, ``t`` linear in
    [0, 1] -- see the derivation above for why this has an exactly-axial tangent at t=0 and
    an exactly-radial tangent at t=1 for any ``shape > 1``.
    """
    if shape <= 1.0:
        raise ValueError(f"shape must be > 1 (endpoint tangent conditions require it); got {shape}")
    if n_points < 2:
        raise ValueError(f"n_points must be >= 2; got {n_points}")
    a = x_te - x_le
    b = r_te - r_le
    if a <= 0.0:
        raise ValueError(f"x_te must be > x_le (the turn must move downstream); got A={a}")
    if b <= 0.0:
        raise ValueError(f"r_te must be > r_le (the turn must expand the radius); got B={b}")

    t = np.linspace(0.0, 1.0, n_points)
    x = x_le + a * t
    # 1 - t**shape lies in [0, 1] for t in [0, 1], so this never hits a negative base --
    # no domain error, no RuntimeWarning, for any shape > 0.
    r = r_te - b * np.power(1.0 - np.power(t, shape), 1.0 / shape)
    return x, r


def _dr_dt(t: float, b: float, shape: float) -> float:
    """Analytic dr/dt of the superellipse turn (see module docstring).

    ``dr/dt = B * t**(shape-1) * (1 - t**shape)**(1/shape - 1)``.

    Returned as the exact IEEE limit at the endpoints (0.0 at t=0, math.inf at t=1) rather
    than evaluated by the formula there, which would either be exactly right anyway (t=0,
    for shape > 1) or hit ``0.0 ** negative`` (t=1), a domain error in plain Python.
    """
    if t <= 0.0:
        return 0.0
    if t >= 1.0:
        return math.inf
    return b * t ** (shape - 1.0) * (1.0 - t**shape) ** (1.0 / shape - 1.0)


def _tangent_slopes(t: float, a: float, b: float, shape: float) -> Tuple[float, float]:
    """(dr/dx, dx/dr) of the superellipse turn at parameter t, as exact IEEE limits.

    t=0 (LE): dr/dt=0    -> (dr/dx, dx/dr) = (0.0, inf)   -- axial tangent.
    t=1 (TE): dr/dt=inf  -> (dr/dx, dx/dr) = (inf, 0.0)   -- radial tangent.

    Used by the test suite to assert the endpoint tangent conditions EXACTLY (``== 0.0``),
    not approximately -- the closed form gives exact values, so a test that only checks
    "close to zero" would be throwing that away.
    """
    dr_dt = _dr_dt(t, b, shape)
    if dr_dt == 0.0:
        return 0.0, math.inf
    if math.isinf(dr_dt):
        return math.inf, 0.0
    return dr_dt / a, a / dr_dt


# --------------------------------------------------------------------------- build_flowpath


def build_flowpath(
    r1h: float,
    r1s: float,
    r2: float,
    b2: float,
    Lz: float,
    r_exit: float,
    inlet_duct_len: float,
    shape: float = 2.0,
    n_points: int = 200,
    lz_refers_to: LzRef = "hub",
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Reconstruct a hub and a shroud meridional curve from a meanline scalar table.

    Parameters
    ----------
    r1h, r1s : inducer hub / shroud radius at the impeller LE, m.
    r2 : impeller exit radius, m.
    b2 : impeller exit passage width, measured ALONG X (docs/PHYSICS-RULES.md rule 4), m.
    Lz : impeller axial length, LE plane -> TE plane, m. Ambiguous between hub and shroud
        (see module docstring) -- which one it means is ``lz_refers_to``.
    r_exit : radius the flow path is carried out to (e.g. a vaneless-diffuser measurement
        plane). The diffuser between r2 and r_exit has PARALLEL walls (constant b = b2) --
        pass ``r_exit = r2`` to omit the diffuser section entirely (used by
        ``eckardt_rotor_o_lm_check`` to isolate the impeller-only meridional length).
    inlet_duct_len : straight axial extension upstream of the LE, m. Pass 0.0 to omit it.
    shape : the superellipse exponent, > 1. THE free reconstruction parameter -- see module
        docstring. 2.0 = ellipse quadrant (default).
    n_points : sample count for EACH of the hub and shroud impeller-turn curves. Straight
        (duct / diffuser) segments always use exactly 2 points (exact for a line).
    lz_refers_to : which curve's LE-to-TE axial extent equals ``Lz``. Default ``"hub"`` --
        see module docstring for the Eckardt-rotor-O justification.

    Returns
    -------
    (hub_xr, shroud_xr) : each an (N, 2) array of (x, r) pairs, m, LE (or inlet-duct start)
        to TE (or diffuser end), suitable for ``write_flowpath_csv`` and, from there,
        ``turbodesign.centrifugal.geometry.MeridionalPath.from_csv``.
    """
    if not (0.0 < r1h < r1s):
        raise ValueError(
            f"require 0 < r1h < r1s (hub below shroud at the LE); got r1h={r1h}, r1s={r1s}"
        )
    if r1s >= r2:
        raise ValueError(
            f"require r1s < r2 (the impeller must expand the radius); got r1s={r1s}, r2={r2}"
        )
    if b2 <= 0.0:
        raise ValueError(f"b2 must be > 0; got {b2}")
    if Lz <= 0.0:
        raise ValueError(f"Lz must be > 0; got {Lz}")
    if r_exit < r2:
        raise ValueError(
            f"r_exit must be >= r2 (can't carry the flow path backward); got r_exit={r_exit} < r2={r2}"
        )
    if inlet_duct_len < 0.0:
        raise ValueError(f"inlet_duct_len must be >= 0; got {inlet_duct_len}")
    if lz_refers_to not in ("hub", "shroud"):
        raise ValueError(f"lz_refers_to must be 'hub' or 'shroud'; got {lz_refers_to!r}")

    x_le = 0.0
    # Sign convention from data/hecc/flowpath_{hub,shroud}.csv (module docstring): the hub
    # TE always sits b2 further downstream (larger x) than the shroud TE, regardless of
    # which curve `Lz` is anchored to.
    if lz_refers_to == "hub":
        x_te_hub = Lz
        x_te_shroud = Lz - b2
        if x_te_shroud <= x_le:
            raise ValueError(
                f"Lz - b2 = {x_te_shroud} <= 0: Lz={Lz} is too small relative to b2={b2} "
                "for lz_refers_to='hub' (the shroud TE would fall upstream of the LE)"
            )
    else:
        x_te_shroud = Lz
        x_te_hub = Lz + b2

    hub_x, hub_r = _superellipse_turn(x_le, r1h, x_te_hub, r2, shape, n_points)
    shroud_x, shroud_r = _superellipse_turn(x_le, r1s, x_te_shroud, r2, shape, n_points)

    if inlet_duct_len > 0.0:
        hub_x = np.concatenate(([x_le - inlet_duct_len], hub_x))
        hub_r = np.concatenate(([r1h], hub_r))
        shroud_x = np.concatenate(([x_le - inlet_duct_len], shroud_x))
        shroud_r = np.concatenate(([r1s], shroud_r))

    if r_exit > r2:
        hub_x = np.concatenate((hub_x, [x_te_hub]))
        hub_r = np.concatenate((hub_r, [r_exit]))
        shroud_x = np.concatenate((shroud_x, [x_te_shroud]))
        shroud_r = np.concatenate((shroud_r, [r_exit]))

    hub_xr = np.column_stack([hub_x, hub_r])
    shroud_xr = np.column_stack([shroud_x, shroud_r])
    return hub_xr, shroud_xr


def write_flowpath_csv(
    name: str,
    hub_xr: npt.NDArray[np.float64],
    shroud_xr: npt.NDArray[np.float64],
    out_dir: str | Path = ".",
) -> Tuple[Path, Path]:
    """Write ``<name>_hub.csv`` / ``<name>_shroud.csv`` (header ``x_m,r_m``, SI metres --
    the NASA HECC Appendix C convention ``MeridionalPath.from_csv`` reads)."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    hub_path = out_path / f"{name}_hub.csv"
    shroud_path = out_path / f"{name}_shroud.csv"
    for path, xr in ((hub_path, hub_xr), (shroud_path, shroud_xr)):
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["x_m", "r_m"])
            writer.writerows(xr.tolist())
    return hub_path, shroud_path


# --------------------------------------------------------------------------- independent check


def meridional_length(hub_xr: npt.NDArray[np.float64], shroud_xr: npt.NDArray[np.float64]) -> float:
    """Arc length of the MEAN line between two meridional curves, start to end, m.

    Hub and shroud curves are, in general, sampled at different point counts and different
    physical arc lengths (this module's own ``build_flowpath`` uses the SAME t-grid for
    both, but a curve read from an arbitrary CSV need not, and this function makes no such
    assumption). They are paired by NORMALISED cumulative arc-length fraction -- point i on
    the hub is matched to the point at the same fraction-of-the-way-along-its-own-curve on
    the shroud, not to the same index and not to the same x or r. The mean line is then
    ``((x_hub+x_shroud)/2, (r_hub+r_shroud)/2)`` sampled on a common fraction grid, and its
    length is the sum of consecutive chord lengths (piecewise-linear approximation -- see
    the discretisation-error discussion in ``tests/unit/test_flowpath_builder.py``).
    """
    if hub_xr.shape[0] < 2 or shroud_xr.shape[0] < 2:
        raise ValueError("each curve needs at least 2 points")

    def _arc_fraction(xr: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        d = np.hypot(np.diff(xr[:, 0]), np.diff(xr[:, 1]))
        s = np.concatenate(([0.0], np.cumsum(d)))
        if s[-1] <= 0.0:
            raise ValueError("curve has zero arc length")
        return s / s[-1]

    f_hub = _arc_fraction(hub_xr)
    f_shroud = _arc_fraction(shroud_xr)

    n = max(hub_xr.shape[0], shroud_xr.shape[0])
    f = np.linspace(0.0, 1.0, n)
    x_hub = np.interp(f, f_hub, hub_xr[:, 0])
    r_hub = np.interp(f, f_hub, hub_xr[:, 1])
    x_shroud = np.interp(f, f_shroud, shroud_xr[:, 0])
    r_shroud = np.interp(f, f_shroud, shroud_xr[:, 1])

    x_mean = 0.5 * (x_hub + x_shroud)
    r_mean = 0.5 * (r_hub + r_shroud)
    return float(np.sum(np.hypot(np.diff(x_mean), np.diff(r_mean))))


# Kovář et al. 2021, *Energies* 14(24) 8545, Table 1 ("Algorithm input parameters"), p.13/22
# -- transcribed and cross-checked in research/notes/07-eckardt-operating-point.md Section 5.
# Diameters in the source converted to radii (r = D/2). UNVERIFIED beyond that transcription
# (the two Eckardt originals are paywalled; Kovář et al. 2021 is the only text-extractable open source
# located for a complete numeric table -- see data/eckardt/README.md).
ECKARDT_ROTOR_O = {
    "r1h": 0.045,  # D1h = 90 mm
    "r1s": 0.140,  # D1s = 280 mm
    "r2": 0.200,  # D2 = 400 mm
    "b2": 0.026,
    "Lz": 0.130,
    "Lm_published": 0.1739,  # Kovář et al. Table 1's own published meridional channel length
}


def eckardt_rotor_o_lm_check(shape: float = 2.0, n_points: int = 2000) -> dict:
    """Compare the RECONSTRUCTED impeller-only meridional length against Kovář et al.'s PUBLISHED
    ``Lm`` for Eckardt rotor O, and report which ``shape`` would reproduce it exactly.

    This is the independent check the module docstring promises: it is not fitting `shape`
    to the published number and then calling that agreement validation (docs/PHYSICS-RULES.md,
    "Failures are reported, not tuned away") -- it reports the DEFAULT reconstruction's
    disagreement, and SEPARATELY reports the shape that would remove it, so both numbers are
    visible. Only the impeller turn is built (``r_exit=r2``, ``inlet_duct_len=0``) because
    Kovář et al.'s ``Lm`` is the meridional CHANNEL length -- LE to TE -- not the inlet duct or the
    downstream diffuser.
    """
    p = ECKARDT_ROTOR_O

    def lm_at(s: float) -> float:
        hub_xr, shroud_xr = build_flowpath(
            p["r1h"],
            p["r1s"],
            p["r2"],
            p["b2"],
            p["Lz"],
            r_exit=p["r2"],
            inlet_duct_len=0.0,
            shape=s,
            n_points=n_points,
            lz_refers_to="hub",
        )
        return meridional_length(hub_xr, shroud_xr)

    lm_default = lm_at(shape)
    disagreement_pct = 100.0 * (lm_default - p["Lm_published"]) / p["Lm_published"]

    # lm_at(shape) is monotonically increasing in shape (shape->1 approaches the straight
    # chord, the shortest possible connector; shape->inf approaches the L-shaped bounding
    # box, the longest) -- bracket generously and let brentq fail loudly if that assumption
    # is ever violated by a future change to _superellipse_turn.
    from scipy.optimize import brentq

    def residual(s: float) -> float:
        return lm_at(s) - p["Lm_published"]

    shape_matching_published_lm = brentq(residual, 1.001, 50.0, xtol=1e-6)

    return {
        "shape_default": shape,
        "Lm_reconstructed_default_m": lm_default,
        "Lm_published_m": p["Lm_published"],
        "disagreement_pct": disagreement_pct,
        "shape_matching_published_Lm": shape_matching_published_lm,
    }


# --------------------------------------------------------------------------- CLI


def _build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct a hub/shroud meridional flow-path CSV pair from a meanline "
            "scalar table (see module docstring for the superellipse derivation and the "
            "Lz / shape reconstruction degrees of freedom)."
        )
    )
    parser.add_argument("--name", help="output basename; writes <name>_hub.csv / <name>_shroud.csv")
    parser.add_argument("--r1h", type=float, help="inducer hub radius, m")
    parser.add_argument("--r1s", type=float, help="inducer shroud radius, m")
    parser.add_argument("--r2", type=float, help="impeller exit radius, m")
    parser.add_argument("--b2", type=float, help="impeller exit passage width (along x), m")
    parser.add_argument("--lz", type=float, help="impeller axial length, LE->TE, m")
    parser.add_argument("--r-exit", type=float, help="flow-path exit radius, m")
    parser.add_argument(
        "--inlet-duct-len", type=float, default=0.0, help="upstream axial duct length, m"
    )
    parser.add_argument(
        "--shape",
        type=float,
        default=2.0,
        help="superellipse exponent, > 1 (default: 2.0, ellipse)",
    )
    parser.add_argument(
        "--n-points", type=int, default=200, help="sample count per impeller-turn curve"
    )
    parser.add_argument(
        "--lz-refers-to",
        choices=("hub", "shroud"),
        default="hub",
        help="which curve Lz is measured along",
    )
    parser.add_argument("--out-dir", default=".", help="output directory")
    parser.add_argument(
        "--report-eckardt-o",
        action="store_true",
        help="also print the Eckardt-rotor-O reconstructed-vs-published Lm check and exit",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_arg_parser()
    args = parser.parse_args(argv)

    if args.report_eckardt_o:
        report = eckardt_rotor_o_lm_check(shape=args.shape, n_points=max(args.n_points, 2000))
        print("Eckardt rotor O -- reconstructed vs. published Lm (Kovář et al. 2021 Table 1)")
        print(f"  shape (default)              : {report['shape_default']}")
        print(f"  Lm reconstructed (default)   : {report['Lm_reconstructed_default_m']:.6f} m")
        print(f"  Lm published                 : {report['Lm_published_m']:.6f} m")
        print(f"  disagreement                 : {report['disagreement_pct']:+.3f} %")
        print(f"  shape matching published Lm  : {report['shape_matching_published_Lm']:.4f}")
        return

    required = ("name", "r1h", "r1s", "r2", "b2", "lz", "r_exit")
    missing = [f"--{name.replace('_', '-')}" for name in required if getattr(args, name) is None]
    if missing:
        parser.error(f"the following arguments are required: {', '.join(missing)}")

    hub_xr, shroud_xr = build_flowpath(
        r1h=args.r1h,
        r1s=args.r1s,
        r2=args.r2,
        b2=args.b2,
        Lz=args.lz,
        r_exit=args.r_exit,
        inlet_duct_len=args.inlet_duct_len,
        shape=args.shape,
        n_points=args.n_points,
        lz_refers_to=args.lz_refers_to,
    )
    hub_path, shroud_path = write_flowpath_csv(args.name, hub_xr, shroud_xr, out_dir=args.out_dir)
    print(f"wrote {hub_path}")
    print(f"wrote {shroud_path}")


if __name__ == "__main__":
    main()
