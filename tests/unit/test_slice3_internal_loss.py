"""Slice 3 -- one internal loss. Root cause D.

Upstream ships 16 centrifugal loss models in loss/compressor/otac.py -- Jansen
clearance, Daily disc friction, Coppage skin friction, Aungier/Oh recirculation,
Johnston mixing, Stanitz vaneless diffuser. Every one declares LossType.Enthalpy.

And compressor_math.py:200-213 handles only Pressure / Polytropic / Entropy:

    if loss_type == LossType.Pressure and callable(loss_fn):
        row.Yp = loss_fn(row, upstream)
    else:
        row.Yp[:] = 0          # <-- Enthalpy lands HERE, silently

So attaching ImpellerVarious() to a rotor produces a LOSS-FREE, ISENTROPIC impeller.
No error. No warning. The entire centrifugal loss suite is dead code.

(The local check_loss_wiring.py pre-commit guard detects this on the live tree today.)

THE DEFINING INVARIANT OF THIS SLICE
------------------------------------
An INTERNAL loss destroys total pressure but does NOT change the work input:

    P02 falls          <- irreversibility
    T02 unchanged      <- no work was added or removed

That is what makes it *internal*. Slice 4 adds PARASITIC losses, which do the exact
opposite. Upstream's single `Yp` cannot express the difference -- which is why its
efficiency curve cannot roll over at any coefficient value.
"""

import pytest

from turbodesign.centrifugal import (
    Air,
    Impeller,
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


def test_internal_loss_destroys_pressure_without_changing_work(path):
    """THE test for this slice.

    Internal loss -> P02 DOWN, T02 UNCHANGED. If T02 moves, the model has confused an
    internal loss with a parasitic one, and efficiency will never behave correctly.
    """
    lossless = _solve(path, losses=None)
    lossy = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))

    te_0 = lossless.stations.impeller_te
    te_1 = lossy.stations.impeller_te

    assert te_1.P0 < te_0.P0, "an internal loss MUST destroy total pressure"
    # NOT exact (rel=1e-9 was wrong -- see test_internal_loss_adds_no_work_term for
    # the actual exact invariant). Internal loss is coupled into the continuity solve:
    # lower P02 -> lower rho2 -> (same mdot) higher Cm2 -> lower Vt2 (beta2b < 0) ->
    # lower Euler work -> T02 genuinely falls, by design. At the HECC design point
    # this is ~0.70 K on a T02 of ~186 K (~0.4%); rel=5e-3 leaves headroom without
    # papering over a real regression in the coupling.
    assert te_1.T0 == pytest.approx(te_0.T0, rel=5e-3), (
        "an internal loss shifts T02 only through the velocity triangle (Cm2/Vt2 "
        "coupling), not directly -- a shift far outside this band means something "
        "else moved, likely a genuine work term leaking in"
    )


def test_internal_loss_lowers_pressure_ratio_and_efficiency(path):
    lossless = _solve(path, losses=None)
    lossy = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))

    assert lossy.PR < lossless.PR
    assert lossy.eta_poly < 1.0
    assert lossy.eta_poly > 0.5, "one skin-friction term should not destroy half the work"


def test_internal_loss_adds_no_work_term(path):
    """An internal loss destroys pressure. It adds NO WORK TERM to the energy equation.

    It DOES shift T02 slightly -- via the velocity triangle, because a lower rho2
    raises Cm2, and Vt2 = sigma*U2 + Cm2*tan(beta2b) with beta2b < 0 then falls.
    That coupling is real (~0.7 K at the HECC design point) and must NOT be
    suppressed. The invariant is on the energy BOOKKEEPING, not on the number:
    """
    lossy = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))
    # no parasitic loss present -> the actual work equals the Euler work exactly
    assert lossy.work_actual == pytest.approx(lossy.work_euler, rel=1e-9)
    assert lossy.losses.parasitic_total == 0.0


def test_entropy_rises_across_an_internal_loss(path):
    """Tier A, 2nd law. Lossless was ds = 0; a real loss must make ds > 0."""
    lossy = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))
    le, te = lossy.stations.impeller_le, lossy.stations.impeller_te
    assert te.s > le.s


# ------------------------------------------------------------- defect D itself


def test_a_configured_loss_model_cannot_silently_contribute_nothing(path):
    """Defect D, directly. Upstream returns a loss-free impeller when a loss model is
    attached whose LossType the solver does not consume -- no error, no warning.

    A loss that is CONFIGURED must CHANGE THE ANSWER, or it must raise. Silence is not
    an option.
    """
    lossless = _solve(path, losses=None)
    lossy = _solve(path, losses=LossSet([ImpellerSkinFrictionJansen()]))
    assert lossy.PR != pytest.approx(lossless.PR, rel=1e-6), (
        "the loss model contributed NOTHING -- this is exactly defect D"
    )


def test_loss_model_declares_its_kind(path):
    """Every loss must self-identify as internal or parasitic. There is no default:
    guessing wrong inverts the physics.
    """
    assert ImpellerSkinFrictionJansen().kind == "internal"


def test_loss_returns_enthalpy_in_joules_per_kg(path):
    """The currency is enthalpy (J/kg), not a pressure-loss coefficient.

    Upstream's OTAC models already return enthalpy -- it simply never consumed them.
    """
    loss = ImpellerSkinFrictionJansen()
    op = _solve(path, losses=LossSet([loss]))
    dh = op.losses.internal["ImpellerSkinFrictionJansen"]
    assert dh > 0.0
    assert dh < 0.5 * op.stations.impeller_te.U**2, "a single loss cannot exceed the Euler work"
