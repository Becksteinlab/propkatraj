"""
Tests that propkatraj can be imported.
"""

import sys

import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF, DCD
import pytest

import propkatraj


def test_single_frame_regression(tmpdir):
    """Single frame propkatraj call compared against same frame written by
    MDA
    """
    u = mda.Universe(PSF, DCD)
    with tmpdir.as_cwd():
        pkas = propkatraj.get_propka(u, start=0, stop=1)

        # There should be one row and 75 columns
        assert len(pkas) == 1
        assert len(pkas.columns) == 75
