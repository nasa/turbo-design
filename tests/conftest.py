"""Shared fixtures.

Validation tests hit published rig data and are slow; they are opt-in:
    uv run pytest                    # unit + oracle only
    uv run pytest -m validation      # the rig matrix
"""

from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA


def pytest_collection_modifyitems(config, items):
    """Skip tests MARKED ``validation`` unless ``-m validation`` is given.

    ⚠ This used to test ``"validation" in item.keywords``, which also matches the *directory
    name*: every test in ``tests/validation/`` was silently skipped whether or not it carried the
    marker. ``tests/validation/test_frozen_coefficients.py`` -- the guard that closes the
    coefficient freeze -- was collected, skipped, and reported green on its first run.

    A guard that does not run is not a guard. Same failure class as the eight "a term that cannot
    fire cannot be caught being wrong" bugs in this project's record (``docs/PHYSICS-RULES.md``),
    one level up: a TEST that cannot fire cannot catch anything being wrong. Key on the marker.
    """
    if config.getoption("-m"):
        return
    skip = pytest.mark.skip(reason="validation: run with `-m validation`")
    for item in items:
        if item.get_closest_marker("validation") is not None:
            item.add_marker(skip)
