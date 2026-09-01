# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import types
sys.path.insert(0, os.path.abspath('_ext'))
import sphinx_rtd_theme

# Ensure project root is on sys.path so turbodesign can be imported when building docs
sys.path.insert(0, os.path.abspath(".."))

# Minimal mocks for optional heavy dependencies used during autodoc import
mock_modules = {
    "pyturbo": types.ModuleType("pyturbo"),
    "pyturbo.helper": types.ModuleType("pyturbo.helper"),
    "pyturbo.aero": types.ModuleType("pyturbo.aero"),
    "pyturbo.aero.airfoil2D": types.ModuleType("pyturbo.aero.airfoil2D"),
    "cantera": types.ModuleType("cantera"),
    "cantera.composite": types.ModuleType("cantera.composite"),
}
for name, module in mock_modules.items():
    sys.modules.setdefault(name, module)

# Provide minimal attributes used in code paths during import
sys.modules["pyturbo.helper"].line2D = lambda *args, **kwargs: None
sys.modules["pyturbo.helper"].convert_to_ndarray = lambda x, *_, **__: x
sys.modules["pyturbo.helper"].xr_to_mprime = lambda *args, **kwargs: None
# NOTE: conf.py is exec'd by Sphinx in a namespace with no "__name__" global, so
# type(name, bases, {}) here would normally build a class with no __module__ at all
# (CPython's type.__new__ infers __module__ from the calling frame's __name__, and
# silently leaves it unset if that lookup fails). Autodoc later stringifies annotations
# like `Optional[Solution]`, which calls repr() on the class and crashes with
# `AttributeError: __module__` for any class missing that attribute. Set __module__
# explicitly on every dynamically-created mock class to avoid that crash.
sys.modules["pyturbo.aero.airfoil2D"].Airfoil2D = type(
    "Airfoil2D", (), {"__module__": "pyturbo.aero.airfoil2D"}
)
mock_solution = type("Solution", (), {"__module__": "cantera"})
sys.modules["cantera"].Solution = mock_solution  # type: ignore[attr-defined]
sys.modules["cantera.composite"].Solution = mock_solution  # type: ignore[attr-defined]

import turbodesign

# -- Project information -----------------------------------------------------

project = 'Turbo Design'
copyright = '2024, Paht Juangphanich'
author = 'Paht Juangphanich <paht.juangphanich@nasa.gov>'

# The full version, including alpha/beta/rc tags.
# Read from pyproject.toml instead of hardcoding, so the docs never drift out of
# sync with the package version again (this was previously hardcoded to '1.0.6'
# while pyproject.toml had already moved on to 1.4.2).
try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib
with open(os.path.join(os.path.dirname(__file__), "..", "pyproject.toml"), "rb") as _f:
    _pyproject = tomllib.load(_f)
version = _pyproject["project"]["version"]
release = version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon'
]
autodoc_mock_imports = ["pyturbo", "cantera"]
autosummary_generate = True
autodoc_class_signature = 'separated'

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = False
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_keyword = True
napoleon_custom_sections = None

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

srclink_project = 'https://github.com/nasa/turbo-design'
srclink_branch = 'main'
srclink_src_path = '/turbodesign'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_theme_options = {
    'collapse_navigation': False,
    'display_version': True,
    'logo_only': True,
    'navigation_depth': 2,
}


html_static_path = ['_static']
rst_context = {'turbo-design': turbodesign}


def setup(app):
    def skip(app, what, name, obj, skip, options):
        members = [
            '__init__',
            '__repr__',
            '__weakref__',
            '__dict__',
            '__module__',
        ]
        return True if name in members else skip

    def rst_jinja_render(app, docname, source):
        src = source[0]
        rendered = app.builder.templates.render_string(src, rst_context)
        source[0] = rendered

    app.connect('autodoc-skip-member', skip)
    app.connect("source-read", rst_jinja_render)
