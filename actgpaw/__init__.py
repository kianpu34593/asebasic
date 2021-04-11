"""
ACTgpaw
A Python package for autonomus convergence test of GPAW. 
"""

# Add imports here
from .utils import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
