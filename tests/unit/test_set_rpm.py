"""Tests for CompressorSpool.set_rpm / TurbineSpool.set_rpm - the bulk RPM push
needed because solve() never re-reads self.rpm on repeat calls."""

import numpy as np
import pytest

from tests.test_shaft_match import _reference_turbine, _solved_compressor


def test_compressor_set_rpm_updates_self_and_every_row():
    spool = _solved_compressor()
    original_rpm = spool.rpm

    spool.set_rpm(original_rpm * 2.0)

    assert spool.rpm == pytest.approx(original_rpm * 2.0)
    for row in spool.rows:
        assert row.rpm == pytest.approx(original_rpm * 2.0)


def test_compressor_blade_speed_reflects_new_rpm_after_resolve():
    spool = _solved_compressor()
    new_rpm = spool.rpm * 1.5

    spool.set_rpm(new_rpm)
    spool.solve()

    rotor = next(r for r in spool.rows if r.row_type.name == "Rotor")
    expected_omega = new_rpm * np.pi / 30.0
    assert np.allclose(rotor.U, expected_omega * rotor.r, rtol=1e-6)


def test_turbine_set_rpm_updates_self_and_every_row():
    turbine = _reference_turbine(rpm=9549.0)
    turbine.set_rpm(5000.0)

    assert turbine.rpm == pytest.approx(5000.0)
    for row in turbine.rows:
        assert row.rpm == pytest.approx(5000.0)
