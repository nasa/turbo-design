"""Slice 4 -- one parasitic loss, and the internal/parasitic SPLIT.

This is the distinction upstream's single `Yp` cannot express, and it is not
bookkeeping -- it changes the thermodynamics.

    INTERNAL  (incidence, skin friction, clearance, mixing, shock):
        destroys P0.  Does NOT change the work input.   T02 unchanged.

    PARASITIC (disc friction, recirculation, leakage):
        ADDS WORK to the fluid.  Produces NO pressure.  T02 RISES, P02 unchanged.

    eta = (dh_Euler - dh_internal - dh_static) / (dh_Euler + dh_parasitic)
                                                   ^^^^^^^^^^^^^^^^^^^^^^
                                     parasitic in the DENOMINATOR (it adds work)

Lumping them into one loss coefficient makes the efficiency ROLLOVER unreproducible
at ANY coefficient value. That is the direct mechanical explanation for the symptom
in the original validation sheet: eta fell MONOTONICALLY with speed (60.3% at 70% N
down to 54.5% at 105% N) instead of peaking near design. The sheet flagged it as
"likely an artifact". It was.

THIS SLICE ALSO CARRIES THE PSI GATE
------------------------------------
Slip alone (Slice 2) gives psi = 0.776 -- the EULER work. NASA measures 0.81.
The residual ~0.035 is parasitic work. It closes HERE, not in Slice 2.
"""

import pytest

from turbodesign.centrifugal import (
    Air,
    Impeller,
    ImpellerDiscFrictionDaily,
    ImpellerSkinFrictionJansen,
    InletState,
    LossSet,
    MeridionalPath,
    Stage,
    WiesnerSlip,
)

R2 = 0.21581
RPM = 21789.0
MDOT = 4.9269
PSI_MEASURED = 0.81  # NASA/CR-2014-218114/REV1, Table A.4a


@pytest.fixture(scope="module")
def path(data_dir):
    return MeridionalPath.from_csv(
        data_dir / "hecc" / "flowpath_hub.csv",
        data_dir / "hecc" / "flowpath_shroud.csv",
    )


def _solve(path, losses):
    stage = Stage(
        path=path,
        impeller=Impeller(
            n_blades=15,
            n_splitters=15,
            # NASA Figure 18 (CR-2014-218114/REV1, PDF p.36), digitized, span-averaged.
            # Supersedes -37.5 (unsourced midpoint of NASA's prose band). See
            # docs/centrifugal/11-blade-angle-discrepancy.md. No inducer blade angle is
            # supplied in this fixture, so le_blade_thickness would have no effect --
            # not added here.
            backsweep_deg=-36.0,
            r_te=R2,
            splitter_le_r=0.0675,
        ),
        diffuser=None,
        slip=WiesnerSlip(),
        losses=losses,
        fluid=Air(),
    )
    return stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))


# ------------------------------------------------------------- the invariant


def test_parasitic_loss_adds_work_without_adding_pressure(path):
    """THE test for this slice -- the exact mirror image of Slice 3.

    Parasitic loss -> T02 RISES (work added), P02 UNCHANGED (no pressure produced).
    """
    base = _solve(path, losses=None)
    para = _solve(path, losses=LossSet([ImpellerDiscFrictionDaily()]))

    te_0, te_1 = base.stations.impeller_te, para.stations.impeller_te

    assert te_1.T0 > te_0.T0, "a parasitic loss MUST add work (raise T02)"
    assert te_1.P0 == pytest.approx(te_0.P0, rel=1e-9), (
        "a parasitic loss MUST NOT produce pressure -- if P02 moved, "
        "this loss is being treated as internal"
    )


def test_parasitic_loss_raises_the_work_factor(path):
    """psi = Cp*dT0/U2^2, and parasitic work raises T02. So psi MUST rise.

    This is the ONLY mechanism by which psi can exceed the Euler value.
    """
    base = _solve(path, losses=None)
    para = _solve(path, losses=LossSet([ImpellerDiscFrictionDaily()]))
    assert para.psi > base.psi


def test_parasitic_loss_lowers_efficiency(path):
    """More work in, same pressure out -> efficiency falls."""
    base = _solve(path, losses=None)
    para = _solve(path, losses=LossSet([ImpellerDiscFrictionDaily()]))
    assert para.eta_poly < base.eta_poly


def test_internal_and_parasitic_are_opposites(path):
    """Side by side, the two kinds must move P0 and T0 in OPPOSITE, disjoint ways.

    If a model cannot express this, its efficiency curve cannot roll over.
    """
    base = _solve(path, losses=None)
    internal = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))
    parasitic = _solve(path, losses=LossSet([ImpellerDiscFrictionDaily()]))

    b = base.stations.impeller_te
    i = internal.stations.impeller_te
    p = parasitic.stations.impeller_te

    # internal: P0 down, T0 approximately fixed (NOT exact -- see
    # test_slice3_internal_loss.py::test_internal_loss_adds_no_work_term for the
    # actual exact invariant: internal loss is coupled into the continuity solve,
    # so a lower rho2 raises Cm2, which shifts Vt2/T02 through the velocity
    # triangle, ~0.7 K at the HECC design point -- not a direct work term)
    assert i.P0 < b.P0 and i.T0 == pytest.approx(b.T0, rel=5e-3)
    # parasitic: T0 up, P0 fixed exactly (parasitic genuinely never touches pressure)
    assert p.T0 > b.T0 and p.P0 == pytest.approx(b.P0, rel=1e-9)


def test_loss_kinds_are_declared_and_disjoint():
    assert ImpellerSkinFrictionJansen().kind == "internal"
    assert ImpellerDiscFrictionDaily().kind == "parasitic"


# ------------------------------------------------------------- THE PSI GATE


def test_parasitic_work_moves_psi_toward_the_measured_value(path):
    """Slice 4 proves the MECHANISM; Slice 5 proves the MAGNITUDE.

    Disc friction is ONE of three parasitic terms. Recirculation and leakage arrive
    with the full Oh set in Slice 5. So psi cannot reach the measured 0.81 here --
    it must only move TOWARD it.
    """
    euler_only = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))
    with_para = _solve(
        path, losses=LossSet([ImpellerSkinFrictionJansen(), ImpellerDiscFrictionDaily()])
    )
    assert with_para.psi > euler_only.psi  # parasitic work RAISES psi
    assert with_para.psi < PSI_MEASURED  # but cannot reach it alone


def test_efficiency_formula_puts_parasitic_in_the_denominator(path):
    """eta = (dh_Euler - dh_internal - dh_static) / (dh_Euler + dh_parasitic)

    Verify the split is actually wired that way, not merely summed.
    """
    op = _solve(
        path,
        losses=LossSet([ImpellerDiscFrictionDaily(), ImpellerSkinFrictionJansen()]),
    )
    assert op.losses.internal_total > 0.0
    assert op.losses.parasitic_total > 0.0
    # parasitic raises the work actually done on the fluid above the Euler work
    assert op.work_actual > op.work_euler
    assert op.work_actual == pytest.approx(op.work_euler + op.losses.parasitic_total, rel=1e-9)
