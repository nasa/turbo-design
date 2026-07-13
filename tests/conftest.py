# tests/conftest.py
"""Load a converged spool from a real example, without editing the example.

Most example scripts wrap their work in a `main()` that returns None. Rather than fork a
copy of their setup into the test suite -- which would let the fixture and the example
drift apart -- we parse the real file and append `return locals()` to main() in the AST.
The file on disk is never modified.

A few examples (e.g. `EEE-HPT/eee_hpt.py`, `radial-turbine/radial_turbine-1D.py`) are
flat top-level scripts with no `main()` at all. For those there is nothing to patch: the
module's own top-level namespace, after exec, already *is* the converged state, so it is
returned directly.

Several examples write their results (JSON, PNG, pickles) to a path built from `__file__`
-- an absolute path inside `examples/` -- so simply running one overwrites files in the
working tree, including any uncommitted edit a contributor has in that directory. Each
example is therefore copied into a scratch directory and the copy is what runs: `__file__`,
`sys.path` and the working directory all point into the copy. The repository is never
written to, so there is nothing to restore afterwards.
"""

import ast
import contextlib
import io
import os
import shutil
import sys
import tempfile
from pathlib import Path

import matplotlib

# A couple of examples call plt.show()/savefig(). Force a non-interactive backend before
# anything else has a chance to import pyplot, so running the suite never pops a GUI
# window or blocks.
matplotlib.use("Agg", force=True)

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = REPO_ROOT / "examples"


def _patch_main_to_return_locals(tree: ast.Module) -> bool:
    """Append `return locals()` to a top-level `main()`. True if there was one."""
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "main":
            node.body.append(
                ast.Return(
                    value=ast.Call(
                        func=ast.Name(id="locals", ctx=ast.Load()), args=[], keywords=[]
                    )
                )
            )
            ast.fix_missing_locations(tree)
            return True
    return False


def _run(path: Path) -> dict:
    """Run one example inside a private copy of its directory; return its namespace.

    `path` points at the example in the repository. The example's whole directory is
    copied to a scratch location and the copy is executed, so the data files the example
    reads are all present and anything it writes lands in the copy.
    """
    sandbox = Path(tempfile.mkdtemp(prefix="turbo-design-example-"))
    try:
        example_dir = sandbox / path.parent.name
        shutil.copytree(path.parent, example_dir)
        script = example_dir / path.name

        tree = ast.parse(script.read_text(), filename=str(script))
        has_main = _patch_main_to_return_locals(tree)

        ns: dict = {"__name__": "__characterization__", "__file__": str(script)}
        # EEE-HPT's script does `from get_ss_ps import ...` -- a bare top-level import of
        # a sibling file. That resolves automatically only when the interpreter is
        # launched as `python eee_hpt.py`, which prepends the script's own directory to
        # sys.path; exec() here does not do that for us.
        sys.path.insert(0, str(example_dir))
        cwd_before = os.getcwd()
        os.chdir(example_dir)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                exec(compile(tree, str(script), "exec"), ns)
                # `__name__` is not "__main__", so an example guarded by
                # `if __name__ == "__main__": main()` has defined main() but not run it.
                return ns["main"]() if has_main else ns
        finally:
            os.chdir(cwd_before)
            sys.path.remove(str(example_dir))
    finally:
        shutil.rmtree(sandbox, ignore_errors=True)


@pytest.fixture(scope="session")
def example_spool():
    """Run an example once per session and return its converged namespace.

    `name` is a path relative to `examples/`, e.g.
    "mattingly-axial-compressor/example9.1.py" or "EEE-HPT/eee_hpt.py".
    """
    cache: dict = {}

    def _load(name: str) -> dict:
        if name not in cache:
            cache[name] = _run(EXAMPLES / name)
        return cache[name]

    return _load
