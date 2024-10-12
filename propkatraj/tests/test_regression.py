"""
Tests that propkatraj can be imported.
"""

import os

import MDAnalysis as mda
import numpy as np
import pandas as pd
import pytest

from MDAnalysisTests.datafiles import PSF, DCD, GRO, XTC
from numpy.testing import assert_almost_equal, assert_equal

import propkatraj
from propkatraj.tests.datafiles import (PSF_FRAME_ZERO_PKA,
                                        PSF_FRAME_NINETY_PKA,
                                        INDEXERR_FRAME14_GRO,
                                        INDEXERR_FRAME14_XTC,
                                        ATTERR_FRAME1_PDB,
                                        ATTERR_FRAME1_XTC,)


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope='module', params=(
    (INDEXERR_FRAME14_GRO, INDEXERR_FRAME14_XTC),
    (ATTERR_FRAME1_PDB, ATTERR_FRAME1_XTC)
))
def u_badframe(request):
    top, traj = request.param
    return mda.Universe(top, traj)


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


@pytest.mark.parametrize('selection', ['ag', 'default'])
def test_single_frame_regression(tmpdir, u, selection):
    """Single frame PropkaTraj regression test"""
    with tmpdir.as_cwd():
        if selection == 'ag':  # test slection by passing an atomgroup
            pkatraj = propkatraj.PropkaTraj(u.select_atoms('protein'))
        else:  # test selection by using the default select="protein"
            pkatraj = propkatraj.PropkaTraj(u)

        pkatraj.run(stop=1)

        # load reference data
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_ZERO_PKA)

        # test residue numbers
        resnums = pkatraj.results.pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test pka values
        assert_almost_equal(pkatraj.results.pkas.values[0], ref_pkas,
                            decimal=2)


def test_multi_frame_regression(tmpdir, u, glu_ref):
    """Multiframe propkatraj call basic regression test

    Note: slow test
    """
    with tmpdir.as_cwd():
        pkatraj = propkatraj.PropkaTraj(u)
        pkatraj.run(step=10)

        # load reference data for resnum and the last frame pka
        ref_resnums, ref_pkas = pka_from_file(PSF_FRAME_NINETY_PKA)

        # test residue numbers
        resnums = pkatraj.results.pkas.columns.to_numpy()
        assert_equal(resnums, ref_resnums)

        # test final frame pka values
        assert_almost_equal(pkatraj.results.pkas.values[-1], ref_pkas, decimal=2)

        # test one data series: glu 162
        assert_almost_equal(pkatraj.results.pkas[162].values, glu_ref, decimal=2)


@pytest.mark.parametrize('start, stop, step', [
    (None, 5, None), (None, None, 25)
])
def test_start_stop_step(tmpdir, u, start, stop, step):
    """Basic test to make sure the dataframe gets populated with the right
    dimensions"""
    with tmpdir.as_cwd():
        propka = propkatraj.PropkaTraj(u)
        propka.run(start=start, stop=stop, step=step)
        pkas = propka.results.pkas

        start, stop, step = u.trajectory.check_slice_indices(start, stop, step)

        nframes = len(range(start, stop, step))

        # There should be nframes rows and 75 columns
        assert len(pkas) == nframes
        assert len(pkas.columns) == 75

        # index should be time
        assert pkas.index.name == 'time'
        times = np.array(range(start, stop, step), dtype=np.float32) + 1
        assert_almost_equal(pkas.index.values, times, decimal=5)


def test_nofilename(tmpdir):
    """See issue #23"""
    # Load a GRO file and then purposely set the filename to None
    u = mda.Universe(GRO)
    u.transfer_to_memory()
    u.trajectory.filename = None

    with tmpdir.as_cwd():
        _ = propkatraj.PropkaTraj(u).run().results.pkas
        assert os.path.isfile('current.pka')


def test_skipframe_error(tmpdir, u_badframe):
    """Test for Issue #10 - Error raised"""
    with tmpdir.as_cwd():
        with pytest.raises(RuntimeError, match="failure on frame"):
            pkatraj = propkatraj.PropkaTraj(u_badframe)
            pkatraj.run()


@pytest.mark.parametrize('top, traj, framenum', [
    (INDEXERR_FRAME14_GRO, INDEXERR_FRAME14_XTC, 14),
    (ATTERR_FRAME1_PDB, ATTERR_FRAME1_XTC, 1)
])
def test_skipframe_pass(tmpdir, caplog, top, traj, framenum):
    """Test for Issue #10 - passes, logging raise"""
    u = mda.Universe(top, traj)
    with tmpdir.as_cwd():
        wmsg = "failure on frame: {0}".format(framenum)
        with pytest.warns(UserWarning, match=wmsg):
            pkatraj = propkatraj.PropkaTraj(u, skip_failure=True)
            pkatraj.run()

        perc = 1 / u.trajectory.n_frames * 100

        msg = ("number of failed frames = 1\n"
               "percentage failure = {0}\n"
               "failed frames: {1}".format(perc, [framenum]))
        assert msg in caplog.records[0].message
        # times should have been modified in-place
        assert len(pkatraj.times) == (len(u.trajectory) - 1)
        assert u.trajectory[framenum].time not in pkatraj.times


def test_nonprotein_residues(tmpdir):
    """Test for Issue #24 - remove once PDBWriter is fixed"""
    u = mda.Universe(GRO, XTC)
    with tmpdir.as_cwd():
        wmsg = "Non protein atoms passed"
        with pytest.warns(UserWarning, match=wmsg):
            _ = propkatraj.PropkaTraj(u, select="not protein")
