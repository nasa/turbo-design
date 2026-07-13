import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_importing_spool_modules_does_not_pull_in_gui_toolkits():
    """The solver must be importable on a headless machine with no Tk.

    turtle (and the tkinter it wraps) has no place in the import graph of a
    numerical library; pulling it in turns a missing display into an
    import-time crash.
    """
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import turbodesign.compressor_spool, turbodesign.turbine_spool; "
            "assert 'turtle' not in sys.modules, 'turtle was imported'; "
            "assert 'tkinter' not in sys.modules, 'tkinter was imported'",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, (
        f"import check failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )
