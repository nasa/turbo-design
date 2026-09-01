"""Tests for the per-row choke-margin diagnostics (flow_math.update_choke_diagnostics)."""

import json

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.flow_math import update_choke_diagnostics
from turbodesign.isentropic import A_As, choke_margin, mass_flow_function


@pytest.mark.parametrize("gamma", [1.3, 1.4, 1.667])
def test_choke_margin_matches_the_a_as_identity(gamma):
    M = np.linspace(0.05, 0.99, 30)
    assert choke_margin(M, gamma) == pytest.approx(1.0 - 1.0 / A_As(M, gamma), rel=1e-9)


def _stator_row(M) -> BladeRow:
    row = BladeRow(row_type=RowType.Stator)
    row.M = np.atleast_1d(np.asarray(M, dtype=float))
    row.gamma = 1.4
    return row


def _rotor_row(M_rel) -> BladeRow:
    row = BladeRow(row_type=RowType.Rotor)
    row.M_rel = np.atleast_1d(np.asarray(M_rel, dtype=float))
    row.M = np.array([0.0])  # deliberately different, to prove the rotor branch ignores it
    row.gamma = 1.4
    return row


def test_stator_uses_absolute_mach():
    row = _stator_row([0.4, 0.5, 0.6])
    update_choke_diagnostics(row)
    expected = mass_flow_function(row.M, row.gamma)
    assert row.mass_flow_function == pytest.approx(expected, rel=1e-9)
    assert row.choke_margin == pytest.approx(choke_margin(row.M, row.gamma), rel=1e-9)
    assert row.choke_margin_min == pytest.approx(float(np.min(row.choke_margin)), rel=1e-9)


def test_rotor_uses_relative_mach_not_absolute():
    row = _rotor_row([0.3, 0.45, 0.55])
    update_choke_diagnostics(row)
    expected = mass_flow_function(row.M_rel, row.gamma)
    assert row.mass_flow_function == pytest.approx(expected, rel=1e-9)
    # If the rotor branch had used row.M (all zero) instead of M_rel, mass_flow_function
    # would be all zero too - assert it is not.
    assert np.all(np.asarray(row.mass_flow_function) > 0)


def test_unsolved_row_is_left_untouched():
    """A row whose Mach hasn't been solved yet (all zero) must not be stamped with
    bogus diagnostics - the caller relies on this to avoid crashing on partial state."""
    row = _stator_row([0.0, 0.0])
    update_choke_diagnostics(row)
    assert np.all(row.mass_flow_function == 0)
    assert np.all(row.choke_margin == 0)
    assert row.choke_margin_min == 0.0


def test_to_dict_json_round_trips():
    row = _stator_row([0.3, 0.5])
    update_choke_diagnostics(row)
    row.percent_hub_shroud = np.array([0.0, 1.0])
    row.x = np.array([0.0, 0.1])
    row.r = np.array([0.2, 0.2])
    d = row.to_dict()

    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.generic):
                return obj.item()
            return super().default(obj)

    decoded = json.loads(json.dumps(d, cls=NumpyEncoder))
    assert decoded["choke_margin_min"] == pytest.approx(row.choke_margin_min, rel=1e-9)
    assert decoded["mass_flow_function"] == pytest.approx(list(row.mass_flow_function), rel=1e-9)
