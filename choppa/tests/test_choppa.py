"""
Unit and regression test for the choppa package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import choppa


def test_choppa_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "choppa" in sys.modules
