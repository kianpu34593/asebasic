"""
Unit and regression test for the actgpaw package.
"""

# Import package, test suite, and other packages as needed
import actgpaw
import pytest
import sys

def test_actgpaw_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "actgpaw" in sys.modules
