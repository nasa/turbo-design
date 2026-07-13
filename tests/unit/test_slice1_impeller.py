"""Slice 1 -- tracer bullet: a rothalpy-consistent isentropic impeller.

The whole diagnosis, proved or refuted, in one slice.

turbodesign's compressor path assumes P0R2,ideal = P0R1 (compressor_math.py:246), which
is true ONLY when U1 == U2 -- i.e. for an axial machine. For HECC that discards a factor
of ~2.78 in relative total pressure and collapses the impeller pressure ratio.

Geometry is NASA's own, from Appendix C of NASA/CR-2014-218114/REV1 (Design-Intent Hot).
Targets are derived from those coordinates, not guessed.
"""

import math

import pytest

from turbodesign.centrifugal import Air, Choked, Impeller, InletState, MeridionalPath, Stage

# NASA HECC, Medic et al. Table 2 (SI)
R2 = 0.21581  # impeller exit radius, m
B2_NASA = 0.01547  # exit blade height, m (0.609 in)
RPM = 21789.0
MDOT = 4.9269  # inlet corrected, kg/s (Table A.4a design point)
BETA2B = -37.5  # backsweep, deg (negative = backswept)
N_BLADES = 15


@pytest.fixture(scope="module")
def path(data_dir):
    return MeridionalPath.from_csv(
        data_dir / "hecc" / "flowpath_hub.csv",
        data_dir / "hecc" / "flowpath_shroud.csv",
    )


@pytest.fixture(scope="module")
def op(path):
    """Isentropic impeller: NO slip, NO losses. The thinnest possible whole machine."""
    stage = Stage(
        path=path,
        impeller=Impeller(n_blades=N_BLADES, backsweep_deg=BETA2B, r_te=R2),
        diffuser=None,
        slip=None,
        losses=None,
        fluid=Air(),
    )
    return stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))


# ---------------------------------------------------------------- the tracer bullet


def test_isentropic_impeller_recovers_the_rothalpy_consistent_pressure_rise(op):
    """THE test. turbodesign returns PR ~1.90 here; correct physics gives ~6.82.

    Same work input either way. The difference is entirely the discarded P0R term.

    Target is 6.817, not the previously-quoted 6.726: the earlier number was computed
    with ``Impeller.blockage=0.05``, a default whose own justification was circular
    (it cited an exit area, 0.020018 m^2, that already contained the 5% being
    "back-solved" -- see ``components.Impeller``'s docstring). Slice 1 is pure
    geometry with no empirical corrections, so ``blockage=0.0`` and PR=6.817.
    """
    assert op.PR == pytest.approx(6.817, rel=0.02)
    assert op.PR > 5.0, "if PR is ~1.9, the P0R2,is = P0R1 bug is still present"


def test_the_discarded_pressure_term_is_recovered(op):
    """P0R2,is / P0R1 must NOT be 1.0. turbodesign assumes exactly 1.0."""
    ratio = op.stations.impeller_te.P0R_is / op.stations.impeller_le.P0R
    assert ratio == pytest.approx(2.778, rel=0.01)
    assert ratio > 2.0


# ---------------------------------------------------------------- geometry guards
# These would have caught the original b2 = 0 bug on day one.


def test_impeller_exit_area_equals_2_pi_r_b(path):
    """area = 2*pi*r*b -- ONE formula, valid at any meridional inclination.

    The original geometry encoded 1.63e-4 m^2 here: 129x too small, because hub and
    shroud were forced to the same x, setting b2 = 0.
    """
    area = path.station_at_radius(R2).area
    assert area == pytest.approx(2 * math.pi * R2 * B2_NASA, rel=0.02)


def test_exit_blade_height_is_measured_along_x(path):
    """At a RADIAL station, b2 spans x -- not r. Hub and shroud share the radius."""
    st = path.station_at_radius(R2)
    assert st.b == pytest.approx(B2_NASA, abs=0.001)
    assert abs(st.r_hub - st.r_shroud) < 0.001, "hub and shroud must share the radius"
    assert abs(st.x_hub - st.x_shroud) == pytest.approx(B2_NASA, abs=0.001)


def test_area_is_positive_when_the_cut_has_decreasing_x(path):
    """The natural radial-exit ordering (shroud at smaller x than hub).

    flow_math.compute_streamline_areas returns a NEGATIVE area for this case.
    """
    st = path.station_at_radius(R2)
    assert st.x_shroud < st.x_hub, "shroud sits upstream of hub at a radial exit"
    assert st.area > 0


def test_meridional_inclination_is_radial_at_the_exit(path):
    """phi = 90 deg at the impeller exit. Upstream fakes this with a 1e6 curvature hack."""
    assert math.degrees(path.station_at_radius(R2).phi) == pytest.approx(90.0, abs=15.0)


# ---------------------------------------------------------------- physics invariants


def test_rothalpy_is_conserved_across_the_impeller(op):
    """I = h + W^2/2 - U^2/2. The -U^2/2 term vanishes only when U is constant."""
    le, te = op.stations.impeller_le, op.stations.impeller_te
    assert te.rothalpy == pytest.approx(le.rothalpy, rel=1e-9)


def test_blade_speed_differs_between_inlet_and_exit(op):
    """U1 != U2 is the whole point. If these are equal the case is degenerate."""
    le, te = op.stations.impeller_le, op.stations.impeller_te
    assert te.U == pytest.approx(492.4, abs=1.0)
    assert te.U / le.U > 2.0


def test_meridional_velocity_decomposition(op):
    """Vm^2 == Vx^2 + Vr^2.

    compressor_math.py:261 uses Vr = W*sin(phi) instead of Vm*sin(phi), so V
    double-counts the tangential component. Invisible at phi=0; catastrophic at phi=90.
    """
    for st in op.stations:
        assert st.Vm**2 == pytest.approx(st.Vx**2 + st.Vr**2, rel=1e-9)


def test_exit_flow_coefficient_is_physical(op):
    """Cm2/U2 ~ 0.15-0.30. With b2 = 0 this would need Cm2 = 11,042 m/s."""
    te = op.stations.impeller_te
    assert 0.15 < te.Vm / te.U < 0.30


def test_no_slip_means_the_work_factor_is_too_high(op):
    """Perfect blade guidance -> psi = 0.872. Physically impossible; NASA measures 0.81.

    Slice 2 adds slip and brings this down. This test pins the BEFORE state so the
    slip cycle has something to move. Target is 0.8719, not the previously-quoted
    0.8645 -- see the PR test above: the old number baked in a circularly-justified
    ``blockage=0.05``, now removed (``blockage=0.0``).
    """
    assert op.psi == pytest.approx(0.8719, rel=0.02)


def test_entropy_does_not_decrease(op):
    """Isentropic here, so it must not RISE either."""
    le, te = op.stations.impeller_le, op.stations.impeller_te
    assert te.s == pytest.approx(le.s, abs=1e-6)


# ---------------------------------------------------------------- fix pins (physics review)


def _stage(path, **impeller_kwargs):
    kwargs = dict(n_blades=N_BLADES, backsweep_deg=BETA2B, r_te=R2)
    kwargs.update(impeller_kwargs)
    return Stage(
        path=path,
        impeller=Impeller(**kwargs),
        diffuser=None,
        slip=None,
        losses=None,
        fluid=Air(),
    )


def test_continuity_solve_does_not_crash_near_choke(path):
    """FIX 1: the old outward-walking bracket TypeErrors here.

    At HECC 21789 rpm, mdot=7.40 kg/s the geometric grid walk (lo=1.0, hi=40.0,
    growth=1.6) steps over BOTH the subsonic and supersonic inlet roots -- they both
    live inside one grid interval -- and walks on into T1 < 0, where
    ``(T1/T01)**n`` for non-integer n raises ``TypeError: '<' not supported between
    'complex' and 'float'``. Bounding the search on [eps, V_star] (the MERIDIONAL
    sonic velocity) finds the subsonic root directly, with no grid to skip over.
    """
    stage = _stage(path)
    op = stage.solve(mdot=7.40, rpm=21789.0, inlet=InletState(P0=101325.0, T0=288.15))
    le = op.stations.impeller_le
    assert le.Vm == pytest.approx(265.34, abs=0.5)
    a1 = Air().speed_of_sound(le.T)
    assert le.Vm / a1 < 1.0, "the returned inlet root must be on the subsonic branch"


def test_root_is_always_on_the_subsonic_branch(op):
    """FIX 2: the [eps, V_star] bound makes the supersonic root unreachable.

    Pin it explicitly, at both stations, for the ordinary design point (not just the
    near-choke case above).
    """
    fluid = Air()
    le, te = op.stations.impeller_le, op.stations.impeller_te
    assert le.Vm / fluid.speed_of_sound(le.T) < 1.0
    assert te.Vm / fluid.speed_of_sound(te.T) < 1.0


def test_excessive_mass_flow_raises_choked_not_a_crash(path):
    """FIX 1: a station that cannot pass the requested mdot at Mm=1 raises ``Choked``.

    Not a ``RuntimeError`` (the old dead-code fallback), not a silent NaN/complex --
    a named, inspectable exception. Slice 10 needs choke detection anyway; this comes
    free from bounding the continuity solve on the sonic meridional velocity.
    """
    stage = _stage(path)
    with pytest.raises(Choked) as excinfo:
        stage.solve(mdot=30.0, rpm=21789.0, inlet=InletState(P0=101325.0, T0=288.15))
    assert excinfo.value.station == "impeller inlet"
    assert excinfo.value.mdot == 30.0
    assert excinfo.value.mdot_choke < 30.0


def test_blockage_defaults_to_zero(path):
    """FIX 3: no silent exit-area fudge. Slice 1 is pure geometry.

    Reintroducing a nonzero default here is exactly the free knob the physics review
    found: blockage 0.00/0.05/0.08 all pass the PR band, so nothing else in this
    suite would catch it coming back.
    """
    assert Impeller(n_blades=N_BLADES, backsweep_deg=BETA2B, r_te=R2).blockage == 0.0


def test_configuring_a_diffuser_is_not_silently_ignored(path):
    """FIX 4: a configured ``diffuser=`` must not vanish with no error and no warning.

    This is the ``LossType.Enthalpy`` failure pattern in embryo (docs/PHYSICS-RULES.md's "single
    most important thing to know" table): a component is wired in, contributes
    nothing, and the caller has no way to know.
    """
    stage = Stage(
        path=path,
        impeller=Impeller(n_blades=N_BLADES, backsweep_deg=BETA2B, r_te=R2),
        diffuser=object(),
        slip=None,
        losses=None,
        fluid=Air(),
    )
    with pytest.raises(NotImplementedError):
        stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))


def test_explicit_x_le_overrides_the_inducer_eye_heuristic(path, op):
    """FIX 6: ``Impeller(x_le=...)`` is a real seam, not a documented-but-unused knob.

    The lossless slice-1 exit state is algebraically independent of the LE station (no
    losses key off it yet), so PR/psi at an explicit ``x_le`` near the inducer eye must
    match the ``inducer_eye()``-heuristic result to machine precision. This becomes a
    load-bearing choice from slice 4 (incidence loss, inducer relative Mach) -- the
    seam has to exist now even though it is inert today.
    """
    default_x = path.inducer_eye().x_hub
    stage = _stage(path, x_le=default_x - 0.01)
    op_explicit = stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))
    assert op_explicit.PR == pytest.approx(op.PR, rel=1e-9)
    assert op_explicit.psi == pytest.approx(op.psi, rel=1e-9)
