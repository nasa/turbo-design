"""overall_entropy_efficiency() must stay <= 1.0 for a near-lossless machine,
unlike overall_polytropic_efficiency()'s single-gamma approximation."""

from tests.test_shaft_match import _reference_turbine, _solved_compressor


def test_compressor_entropy_efficiency_does_not_exceed_one():
    spool = _solved_compressor()
    assert 0.0 < spool.overall_entropy_efficiency() <= 1.0 + 1e-6


def test_turbine_entropy_efficiency_does_not_exceed_one():
    turbine = _reference_turbine(rpm=9549.0)
    turbine.solve()
    assert 0.0 < turbine.overall_entropy_efficiency() <= 1.0 + 1e-6
