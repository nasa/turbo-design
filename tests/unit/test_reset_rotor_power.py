"""Tests for flow_math.reset_rotor_power - the state-reset primitive needed before
re-solving an already-solved spool at a new operating point."""

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.flow_math import reset_rotor_power


def _row(row_type, power=100.0, power_mean=90.0):
    row = BladeRow(row_type=row_type)
    row.power = power
    row.power_mean = power_mean
    return row


def test_zeroes_rotor_power_only():
    rotor = _row(RowType.Rotor)
    stator = _row(RowType.Stator)
    igv = _row(RowType.IGV)
    inlet = _row(RowType.Inlet)

    reset_rotor_power([rotor, stator, igv, inlet])

    assert rotor.power == 0.0
    assert rotor.power_mean == 0.0
    assert stator.power == 100.0
    assert stator.power_mean == 90.0
    assert igv.power == 100.0
    assert inlet.power == 100.0


def test_empty_list_is_a_no_op():
    reset_rotor_power([])  # must not raise
