"""
BASIC
A Python package for Bulk Adsorption Surface energy calculation with automatIc Convergence
"""

# Add imports here
from .utils import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
