"""Regression test for the stator P0 staleness bug: a lossless compressor stator
must conserve total pressure exactly, not silently drift from a stale seed."""

import pytest

from tests.test_shaft_match import _solved_compressor


def test_lossless_stator_conserves_total_pressure():
    spool = _solved_compressor()
    rows = spool._all_rows()
    rotor, stator = rows[1], rows[2]

    assert stator.P0 == pytest.approx(rotor.P0, rel=1e-6)


def test_lossless_stator_has_no_spurious_negative_entropy_rise():
    spool = _solved_compressor()
    stator = spool._all_rows()[2]

    assert stator.entropy_rise >= -1e-9
