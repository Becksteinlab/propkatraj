# propkatrj --- https://github.com/Becksteinlab/propkatraj
# Copyright (c) 2013-2017 David Dotson, Ricky Sexton, Armin Zjajo, Oliver Beckstein
# Released under the GNU General Public License v3+

from __future__ import print_function

from six import string_types
import os
import tempfile
from six import StringIO, raise_from

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

import warnings
import logging


class PropkaTraj(AnalysisBase):
    """Add some docstring here"""
    def __init__(self, universe, sel='protein', start=None, stop=None,
                 skip_failure=False, **kwargs):
        """Init routine for the PropkaTraj class

        TODO: improve this

        Parameters
        ----------
        universe : :class:`MDAnalysis.Universe` or
                   :class:`MDAnalysis.AtomGroup`
            MDAnalysis Universe or AtomGroup to obtain pKas for.
        sel: str, array_like
            Selection string to use for selecting atoms to use from given
            ``universe``. Can also be a numpy array or list of atom indices.
            [`protein`]
        skip_failure : bool
            If set to ``True``, skip frames where PROPKA fails. If ``False``
            raise an exception. [`False`]

        Results
        -------
        pkas : :class:`pandas.DataFrame`
            DataFrame giving estimated pKa value for each residue for each
            trajectory frame. Residue numbers are given as column labels,
            times as row labels.
        """
        if isinstance(sel, string_types):
            self.ag = universe.select_atoms(sel)
        elif isinstance(sel, (list, np.ndarray)):
            self.ag = universe.atoms[sel]
        else:
            self.ag = universe.atoms

        # probably needs to be a temporary directory instead
        self.tmpfile = os.path.join(os.path.dirname(self.ag.universe.filename),
                                    'current.pdb')
        self.skip_failure = skip_failure
        super(PropkaTraj, self).__init__(self.ag.universe.trajectory, **kwargs)

    def _prepare(self):
        """Creates necessary containers to store pka values"""
        self._pkas = []
        self.num_failed_frames = 0
        self.failed_frames_log = []
        self.failed_times = []
        self._columns = None

    def _single_frame(self):
        pstream = mda.lib.util.NamedStream(StringIO(), self.tmpfile)
        self.ag.write(pstream)
        # reset stream for reading
        pstream.reset()

        try:
            mol = pk.single(pstream, optargs=['--quiet'])
        except (IndexError, AttributeError) as err:
            errmsg = "failure on frame: {0}".format(self._ts.frame)
            if not self.skip_failure:
                raise_from(RuntimeError(errmsg), None)
            else:
                warnings.warn(errmsg)
                self.num_failed_frames += 1
                self.failed_frames_log.append(self._ts.frame)
                self.failed_times.append(self._ts.time)
        else:
            confname = mol.conformation_names[0]
            conformation = mol.conformations[confname]
            groups = conformation.get_titratable_groups()

            # extract pka estimates from each residue
            self._pkas.append([g.pka_value for g in groups])
            if self._columns is None:
                self._columns=[g.atom.resNumb for g in groups]
        finally:
            # deallocate stream
            pstream.close(force=True)

    def _conclude(self):
        # to allow for popping times later, convert self.times to list
        times = self.times.tolist()

        # Ouput failed frames
        if self.num_failed_frames > 0:
            perc_failure = (self.num_failed_frames / self.n_frames) * 100
            # if frames have failed we need to ammend times accordingly
            for i in self.failed_times:
                times.remove(i)

            wmsg = ("number of failed frames = {0}\n"
                    "percentage failure = {1}\n"
                    "failed frames: {2}".format(self.num_failed_frames,
                                                perc_failure,
                                                self.failed_frames_log))
            warnings.warn(wmsg)

        self.pkas = pd.DataFrame(self._pkas, index=pd.Float64Index(times,
                                 name='time'), columns=self._columns)


def get_propka(universe, sel='protein', start=None, stop=None, step=None, skip_failure=False):
    """Get and store pKas for titrateable residues near the binding site.

    Parameters
    ----------
    universe : :class:`MDAnalysis.Universe`
        Universe to obtain pKas for.
    sel : str, array_like
        Selection string to use for selecting atoms to use from given
        ``universe``. Can also be a numpy array or list of atom indices to use.
    start : int
        Frame of trajectory to start from. `None` means start from beginning.
    stop : int
        Frame of trajectory to end at. `None` means end at trajectory end.
    step : int
        Step by which to iterate through trajectory frames. propka is slow,
        so set according to how finely you need resulting timeseries.
    skip_failure : bool
        If set to ``True``, skip frames where PROPKA fails. If ``False``
        raise an exception. The default is ``False``.
        Log file (at level warning) contains information on failed frames.

    Results
    -------
    pkas : :class:`pandas.DataFrame`
        DataFrame giving estimated pKa value for each residue for each
        trajectory frame. Residue numbers are given as column labels, times as
        row labels.

    """
   
    # need AtomGroup to write out for propka
    if isinstance(sel, string_types):
        atomsel = universe.select_atoms(sel)
    elif isinstance(sel, (list, np.ndarray)):
        atomsel = universe.atoms[sel]

    # "filename" for our stream
    # use same name so that propka overwrites
    newname = os.path.join(os.path.dirname(universe.filename), 'current.pdb')

    # progress logging output (because this is slow...)
    pm = mda.lib.log.ProgressMeter(universe.trajectory.n_frames,
                                   format="{step:5d}/{numsteps} t={time:12.3f} ps  "
                                   "[{percentage:5.1f}%]",
                                   interval=1)

    times = []
    pkas = []
    failed_frames = 0
    failed_frames_log = []
    for ts in universe.trajectory[start:stop:step]:
        pm.echo(ts.frame, time=ts.time)

        # we create a named stream to write the atoms of interest into
        pstream = mda.lib.util.NamedStream(StringIO(), newname)
        atomsel.write(pstream)

        pstream.reset()         # reset for reading

        # we feed the stream to propka, and it reads it as if it were a file on
        # disk
        try:
               mol = pk.single(pstream, optargs=['--quiet'])
        except (IndexError, AttributeError) as err:
                     #https://github.com/Becksteinlab/propkatraj/issues/13
                     #https://github.com/Becksteinlab/propkatraj/issues/10
               if not skip_failure:          
                   raise
               else:
                   err_msg = "{0} (failure {2}): failing frame {1}".format(
                         universe.trajectory.filename, ts.frame, failed_frames)
                   failed_frames += 1
                   failed_frames_log.append(ts.frame)
                   logging.warning(err_msg)
                   continue
        finally:
               pstream.close(force=True)  # deallocate

        # parse propka data structures to get out what we actually want
        confname = mol.conformation_names[0]
        conformation = mol.conformations[confname]
        groups = conformation.get_titratable_groups()

        # extract pka estimates from each residue
        pkas.append([g.pka_value for g in groups])

        # record time
        times.append(ts.time)

    if failed_frames_log:
       logging.warning('number of failed frames = {0}'.format(failed_frames))
       logging.warning('percent failure = {0:.3f}%'.format(float(failed_frames)/len(universe.trajectory)*100))
       logging.warning('failed frames: %r', failed_frames_log)
   
    # a `pandas.DataFrame` is a good data structure for this data
    df = pd.DataFrame(pkas, index=pd.Float64Index(times, name='time'),
                      columns=[g.atom.resNumb for g in groups])

    return df
