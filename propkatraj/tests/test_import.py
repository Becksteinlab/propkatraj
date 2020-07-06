"""
Tests that propkatraj can be imported.
"""

import sys
import pytest
import propkatraj


def test_import_basic():
    assert "propkatraj" in sys.modules
