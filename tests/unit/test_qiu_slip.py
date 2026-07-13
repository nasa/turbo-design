"""Slice S6 -- Qiu's real dbeta/dm. Delete the proxy. Do NOT adopt Qiu.

docs/centrifugal/17-tdd-plan.md S6; docs/centrifugal/15-experiments-prereg.md E5;
docs/centrifugal/13-literature-review.md S1.3; data/coefficients.md S5.

THE DEFECT (before this slice). ``QiuSlip.sigma``'s ``dsigma_turn`` substituted
``tan(beta2b)`` for Qiu 2011 Eq. (10b)'s real term ``s2*(dbeta/dm)_2/cos(beta2b)``.
That proxy is driven by the (positive) MAGNITUDE of backsweep, so it is blind to the
one thing Eq. (10b) is actually about -- the SIGN and SIZE of the blade's turning rate
near the exit -- and for HECC's S-shaped blade (Fig 18: |beta| has a mid-chord minimum
and RISES toward the TE) it gets the sign backward: Qiu S3.2 states plainly that a
blade angle that DECREASES toward the exit (signed dbeta/dm < 0, this codebase's
beta2b_deg convention -- negative for backsweep) makes sigma INCREASE with the exit
flow coefficient phi2. HECC's real geometry has signed (dbeta/dm)_2 < 0 (measured,
extract_hecc_blade_angles.py); the proxy inverts the trend regardless.

THE FIX. ``dsigma_turn = F*s2*phi2*(dbeta/dm)_2 / (4*cos(beta2b))``, ``s2 = 2*pi*r2/Z``
(Qiu's own Nomenclature: "s = pitch at the blade exit; s = 2*pi*R2/Z"), with
``dbeta/dm`` and ``r2`` threaded through the ``SlipModel`` protocol's ``**kw`` seam.
The proxy is DELETED, not patched: absent ``dbeta/dm``, ``dsigma_turn = 0`` -- Qiu's
own treatment of an unavailable gradient term (S3.3, for ``dsigma_passage``: "it will
always be assumed to be zero in all of our validation studies").

WHY THIS SLICE MUST BE OUTPUT-INVARIANT. QiuSlip is NOT the default (WiesnerSlip is)
and is not used by any driver or fixture. Every shipped number must stay bit-for-bit
identical -- ``test_shipped_outputs_are_bit_identical`` below, and
``tests/unit/test_ledger_truth.py``, are the guards.
"""

from __future__ import annotations

import inspect
import math
import sys
import warnings
from pathlib import Path

import pytest

from turbodesign.centrifugal import InletState, QiuSlip, WiesnerSlip

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from hecc_stage import BACKSWEEP, IMPELLER, MDOT_DESIGN, RPM, P01, T01, build  # noqa: E402
from test_ledger_truth import (  # noqa: E402
    PSI_AFTER_S5,
    STAGE_ETA_POLY_REALGAS_AFTER_S5,
    STAGE_PR_AFTER_S5,
)

# Real HECC geometry, IMPORTED (not duplicated) from the driver -- docs/centrifugal/
# 17-tdd-plan.md's own rule (b): never duplicate configuration into a test.
Z_EXIT = (
    IMPELLER["n_blades"] + IMPELLER["n_splitters"]
)  # 30 -- BOTH blade rows reach the TE
R2 = IMPELLER["r_te"]
BETA2B_DEG = BACKSWEEP  # -36.0, NASA Fig 18


# --------------------------------------------------------------------- the oracle


def test_sigma_rises_with_flow_coefficient_when_dbeta_dm_is_negative():
    """THE ORACLE, and it is INDEPENDENT of HECC's own measured dbeta/dm.

    Qiu 2011 S3.2, verbatim: "if the blade angle decreases toward the impeller exit,
    d(beta)/dm < 0, the turning term reduces the blade loading and, therefore, causes
    the slip factor to increase with the exit flow coefficient." Eckardt rotor A -- the
    same S-shaped blade family as HECC -- has a MEASURED slip factor that RISES with
    phi2; Qiu 2011 notes Wiesner's and Stodola's models both failed badly for that case
    because neither responds to phi2 at all.

    THE PROXY INVERTED THIS (see module docstring): with the proxy, sigma FALLS as
    phi2 rises for any backswept blade, because tan(abs(beta2b)) > 0 regardless of the
    blade's actual turning direction. This test FAILS on the pre-slice-S6 code.
    """
    q = QiuSlip()
    dbeta_dm = -3.5  # signed, rad/m -- Qiu S3.2's own sign convention; magnitude is
    # immaterial to this test, only that it is NEGATIVE (a decreasing blade angle
    # toward the exit, HECC's own S-shaped-blade case, Fig 18)
    sigma_low = q.sigma(
        beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.10, dbeta_dm=dbeta_dm, r2=R2
    )
    sigma_high = q.sigma(
        beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.35, dbeta_dm=dbeta_dm, r2=R2
    )
    assert sigma_high > sigma_low, (
        "with signed dbeta/dm < 0, sigma MUST rise with the exit flow coefficient "
        "(Qiu 2011 S3.2) -- if it falls, the turning term's sign is wrong"
    )


def test_sigma_at_hecc_with_real_geometry():
    """The design-point oracle, docs/centrifugal/15-experiments-prereg.md E5's own
    "Qiu FULL, real dbeta/dm" row: Z=30 (BOTH blade rows reach the exit -- see the
    class docstring's point 4), dbeta/dm=-3.5 rad/m (the pre-registered figure), on
    HECC's own backsweep and exit radius. sigma = 0.9327 +- 0.005.

    This is a FORMULA-correctness check against a hand-computable oracle, not a claim
    about the exact re-derived dbeta/dm (which the code now carries via
    Impeller.dbeta_dm_te, and which turns out to differ somewhat from -3.5 -- see
    data/coefficients.md S5 and Impeller.dbeta_dm_te's own docstring for the honest
    re-derivation and its window sensitivity).
    """
    sigma = QiuSlip().sigma(
        beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.18, dbeta_dm=-3.5, r2=R2
    )
    assert sigma == pytest.approx(0.9327, abs=0.005)


# --------------------------------------------------------------------- the proxy is gone


def test_the_wrong_signed_proxy_is_unreachable():
    """No code path returns the deleted ``tan(beta2b)`` proxy form."""
    src = inspect.getsource(QiuSlip.sigma)
    assert "math.tan(beta2b)" not in src, (
        "the deleted proxy expression must not survive"
    )
    assert "phi2 * math.tan" not in src, "the deleted proxy expression must not survive"


def test_absent_gradient_makes_sigma_independent_of_flow_coefficient():
    """Behavioural companion to the source-inspection check above: with no ``dbeta_dm``
    supplied at all (the caller does not have the geometry), sigma must not vary with
    phi2 -- the OLD proxy's whole reason for existing was to manufacture exactly this
    variation from nothing. Absent the real gradient, there is nothing left to vary.
    """
    q = QiuSlip()
    s_low = q.sigma(beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.10)
    s_high = q.sigma(beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.35)
    assert s_low == pytest.approx(s_high, rel=1e-12)


# --------------------------------------------------------------------- Qiu's own reductions


def test_no_gradient_reduces_to_the_radial_term():
    """Qiu's own treatment of an unavailable gradient term (S3.3, stated for
    dsigma_passage, applied here to dsigma_turn absent dbeta/dm): dropping it means
    sigma = 1 - dsigma_radial exactly, whatever phi2 is.
    """
    beta2b = math.radians(abs(BETA2B_DEG))
    gamma2 = math.radians(90.0)
    F = 1.0 - 2.0 * math.sin(math.pi / Z_EXIT) * math.sin(
        math.pi / Z_EXIT + beta2b
    ) * math.cos(beta2b) * math.sin(gamma2)
    d_radial = F * math.pi * math.cos(beta2b) * math.sin(gamma2) / Z_EXIT
    expected = 1.0 - d_radial

    sigma = QiuSlip().sigma(beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.25)
    assert sigma == pytest.approx(expected, rel=1e-12)


def test_stodola_reduction_is_preserved():
    """gamma2=90, F->1 (Z->inf), dbeta/dm=0 -> sigma -> 1 - pi*cos(beta2b)/Z.

    Unchanged by this slice (only dsigma_turn's closure changed) -- a regression guard
    that S6 did not disturb the existing radial-term / Stodola sanity check.
    """
    beta2b_deg = -20.0
    Z = 1.0e6  # F -> 1 in this limit; sin(pi/Z) -> 0
    sigma = QiuSlip().sigma(
        beta2b_deg=beta2b_deg, Z=Z, Cm2_over_U2=0.3, dbeta_dm=0.0, r2=0.2
    )
    expected = 1.0 - math.pi * math.cos(math.radians(abs(beta2b_deg))) / Z
    assert sigma == pytest.approx(expected, rel=1e-6)


# --------------------------------------------------------------------- output invariance


def test_shipped_outputs_are_bit_identical():
    """THE constraint on this whole slice. WiesnerSlip is the default in every driver
    and fixture; QiuSlip is unused. Wiring dbeta_dm/r2 through Stage.solve's **kw seam
    must not move a single shipped digit.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        op = build(BACKSWEEP).solve(
            mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )

    assert op.psi == pytest.approx(PSI_AFTER_S5, rel=1e-12)
    assert op.stage_PR == pytest.approx(STAGE_PR_AFTER_S5, rel=1e-12)
    assert op.stage_eta_poly_realgas == pytest.approx(
        STAGE_ETA_POLY_REALGAS_AFTER_S5, rel=1e-12
    )


def test_wiesner_ignores_the_new_kwargs_entirely():
    """WiesnerSlip's ``**kw`` seam must swallow ``dbeta_dm``/``r2`` with ZERO effect --
    the mechanism ``Stage.solve`` relies on for output invariance, exercised directly.
    """
    w = WiesnerSlip()
    plain = w.sigma(beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.2)
    with_extra = w.sigma(
        beta2b_deg=BETA2B_DEG, Z=Z_EXIT, Cm2_over_U2=0.2, dbeta_dm=-99.0, r2=0.5
    )
    assert plain == with_extra
