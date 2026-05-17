import os
import sys
from pathlib import Path
import warnings

# Check if we're building on ReadTheDocs
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if on_rtd:
    # On ReadTheDocs, the package is installed in the environment
    # No need to modify sys.path
    pass
else:
    # Local development - add path to package
    sys.path.insert(0, os.path.abspath('..'))
    sys.path.insert(0, os.path.abspath('../../'))

# -- Project information -----------------------------------------------------
project = 'FRB Repository'
copyright = '2025'
author = 'The FRB Community'

# The full version, including alpha/beta/rc tags
release = '2.3.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary', 
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
suppress_warnings = ['ref.python']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# -- Extension configuration -------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_type_aliases = None
napoleon_custom_sections = [
    ('Args', 'params_style'),
]

# Keep docs build focused on documentation issues rather than runtime noise
warnings.filterwarnings(
    'ignore',
    message='Please define the variable EAZYDIR in your environment pointing to the EAZY folder.',
    category=UserWarning,
)
warnings.filterwarnings('ignore', category=SyntaxWarning, message='invalid escape sequence.*')
warnings.filterwarnings('ignore', message='more than one target found for cross-reference .*')

try:
    from astropy.utils.exceptions import AstropyDeprecationWarning
except Exception:  # pragma: no cover
    AstropyDeprecationWarning = Warning
warnings.filterwarnings('ignore', category=AstropyDeprecationWarning)

# Intersphinx configuration
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'xarray': ('https://xarray.pydata.org/en/stable/', None),
}

# -- Options for autodoc ----------------------------------------------------

# Mock imports for optional dependencies so autodoc can still document modules
# that depend on packages not installed in the build environment.
autodoc_mock_imports = [
    'pcigale',
    'ppxf',
    'pymc3',
    'threedhst',
    'pathos',
    'hmf_emulator',
    'pyregion',
    'spectral_cube',
    'specdb',
    'dl',           # datalab-client
    'pyvo',
    'photutils',
    'dust_extinction',
    'scikit_image',
    'skimage',
    'astroquery',
    'pysftp',
    'asymmetric_kde',
    'theano',
    'IPython',
    'ne2001',
    'pygedm',
    'pdf_fns',   # local module in frb/dm_kde imported as bare name
    'sklearn',
]

# Provide stub env vars so modules that check them at import time don't raise
import os as _os
_os.environ.setdefault('FRB_GDB', '/tmp/stub_gdb')
_os.environ.setdefault('FRB_DATA', '/tmp/stub_data')

autodoc_default_options = {
    'show-inheritance': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'exclude-members': '__weakref__'
}