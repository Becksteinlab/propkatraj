"""
Tests that propkatraj can be imported.
"""
from __future__ import print_function

import sys

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from MDAnalysisTests.datafiles import PSF, DCD
from numpy.testing import assert_almost_equal, assert_equal

import propkatraj
from propkatraj.tests.datafiles import PSF_FRAME_ZERO_PKA, PSF_FRAME_NINETY_PKA


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope='module')
def glu_ref():
    arr = np.array([2.53, 2.61, 2.84, 2.57, 2.76, 2.51, 2.68, 2.95,
                    3.18, 3.34])
    return arr


def pka_from_file(filename):
    columns = ['resname', 'resnum', 'chain', 'pka', 'model-pka']
    df = pd.read_csv(filename, names=columns, skiprows=1,
                     delim_whitespace=True)
    df.sort_values(by=['resnum'], inplace=True)
    resnums = df['resnum'].to_numpy()
    pkas = df['pka'].to_numpy()
    return resnums, pkas


@pytest.mark.parametrize('selection', ['protein', 'array'])
def test_single_frame_regression(tmpdir, u, selection):
    """Single frame propkatraj call compared against same frame written by
    MDA
    """
    with tmpdir.as_cwd():
        if selection == 'array':
            selection = u.select_atoms('protein').ix

        pkas = propkatraj.get_propka(u, sel=selection, stop=1)

        # load reference data
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_ZERO_PKA)

        # test residue numbers
        resnums = pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test pka values
        assert_almost_equal(pkas.values[0], ref_pkas, decimal=2)


@pytest.mark.parametrize('selection', ['protein', 'ag', 'array'])
def test_single_frame_regression_analysisbase(tmpdir, u, selection):
    """Single frame PropkaTraj regression test"""
    with tmpdir.as_cwd():
        if selection == 'ag':
            pkatraj = propkatraj.PropkaTraj(u.select_atoms('protein'),
                                            sel=None)
        elif selection == 'array':
            selection = u.select_atoms('protein').ix
            pkatraj = propkatraj.PropkaTraj(u, sel=selection)
        else:
            pkatraj = propkatraj.PropkaTraj(u)

        pkatraj.run(stop=1)

        # load reference data
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_ZERO_PKA)

        # test residue numbers
        resnums = pkatraj.pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test pka values
        assert_almost_equal(pkatraj.pkas.values[0], ref_pkas, decimal=2)


def test_multi_frame_regression(tmpdir, u, glu_ref):
    """Multiframe propkatraj call basic regression test

    Note: slow test
    """
    with tmpdir.as_cwd():
        pkas = propkatraj.get_propka(u, step=10)

        # load reference data for resnum and last frame pka
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_NINETY_PKA)

        # test residue numbers
        resnums = pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test final frame pka values
        assert_almost_equal(pkas.values[-1], ref_pkas, decimal=2)

        # test one data series: glu 162
        assert_almost_equal(pkas[162].values, glu_ref, decimal=2)


def test_multi_frame_regression_analysisbase(tmpdir, u, glu_ref):
    """Multiframe propkatraj call basic regression test

    Note: slow test
    """
    with tmpdir.as_cwd():
        pkatraj = propkatraj.PropkaTraj(u)
        pkatraj.run(step=10)

        # load reference data for resnum and the last frame pka
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_NINETY_PKA)

        # test residue numbers
        resnums = pkatraj.pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test final frame pka values
        assert_almost_equal(pkatraj.pkas.values[-1], ref_pkas, decimal=2)

        # test one data series: glu 162
        assert_almost_equal(pkatraj.pkas[162].values, glu_ref, decimal=2)


@pytest.mark.parametrize('start, stop, step', [
    (None, 5, None), (None, None, 25)
])
def test_start_stop_step(tmpdir, u, start, stop, step):
    """Basic test to make sure the dataframe gets populated with the right
    dimensions"""
    with tmpdir.as_cwd():
        pkas = propkatraj.get_propka(u, start=start, stop=stop, step=step)

        start, stop, step = u.trajectory.check_slice_indices(start, stop, step)

        nframes = len(range(start, stop, step))

        # There should be nframes rows and 75 columns
        assert len(pkas) == nframes
        assert len(pkas.columns) == 75

        # index should be time
        assert pkas.index.name == 'time'
        times = np.array(range(start, stop, step), dtype=np.float32) + 1
        assert_almost_equal(pkas.index.values, times, decimal=5)


@pytest.mark.parametrize('start, stop, step', [
    (None, 5, None), (None, None, 25)
])
def test_start_stop_step_analysisbase(tmpdir, u, start, stop, step):
    """Basic test to make sure the dataframe gets populated with the right
    dimensions"""
    with tmpdir.as_cwd():
        pkatraj = propkatraj.PropkaTraj(u)
        pkatraj.run(start=start, stop=stop, step=step)

        start, stop, step = u.trajectory.check_slice_indices(start, stop, step)
        nframes = len(range(start, stop, step))

        # There should be nframes rows and 75 columns
        assert len(pkatraj.pkas) == nframes
        assert len(pkatraj.pkas.columns) == 75

        # index should be time
        assert pkatraj.pkas.index.name == 'time'
        times = np.array(range(start, stop, step), dtype=np.float32) + 1
        assert_almost_equal(pkatraj.pkas.index.values, times, decimal=5)
