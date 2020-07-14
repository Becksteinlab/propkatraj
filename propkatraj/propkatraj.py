# propkatrj --- https://github.com/Becksteinlab/propkatraj
# Copyright (c) 2013-2020 David Dotson, Irfan Alibay, Ricky Sexton,
#                         Armin Zjajo, Shujie Fan, Oliver Beckstein
# Released under the GNU General Public License v3+

from __future__ import print_function, division

from six import string_types, StringIO, raise_from
import os

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.util import deprecate

import warnings
import logging


class PropkaTraj(AnalysisBase):
    """Per residue pKa analysis of a trajectory.

    Runs :program:`propka` on the titrateable residues of the selected
    AtomGroup on each frame in the trajectory. Run the analysis with
    :meth:`PropkaTraj.run`, and pKa values will be stored in a
    :class:`pandas.DataFrame` named :attr:`PropkaTraj.pkas`.

    Parameters
    ----------
    atomgroup : :class:`MDAnalysis.Universe` or :class:`MDAnalysis.AtomGroup`
        Group of atoms containing the residues for pKa analysis. Please note
        that :class:`MDAnalysis.UpdatingAtomGroup` are not supported and will
        be automatically converted to :class:`MDAnalysis.AtomGroup`.
    select : str
        Selection string to use for selecting a subsection of atoms to use
        from the input ``atomgroup``. Note: passing non-protein residues to
        :program:`propka` may lead to incorrect results (see notes).
        [`protein`]
    skip_failure : bool
        If set to ``True``, skip frames where :program:`propka` fails. A list
        of failed frames is made available in
        :attr:`PropkaTraj.failed_frames_log`. If ``False`` raise a
        RuntimeError exception on those frames. [`False`]


    Notes
    -----
    Currently only the default behaviour supplemented with the `--quiet` flag
    of :program:`propka` is used.

    Temporary :program:`propka` files are written in the current working
    directory. This will leave a ``current.pka`` and ``current.propka_input``
    file. These are the temporary files for the final frame and can be removed
    safely.

    Current known issues:

    1. Due to the current behaviour of the MDAnalysis PDBWriter, non-protein
       atoms are written to PDBs using `ATOM` records instead of `HETATM`.
       This is likely to lead to undefined behaviour in :program:`propka`,
       which will likely expect `HETATM` inputs. We recommend users to only
       pass protein atoms for now. See the following issue for more details:
       https://github.com/Becksteinlab/propkatraj/issues/24


    Examples
    --------
    A common use case is to calculate the pKa values of titrateable residues in
    a protein across a trajectory.

    To do this you first must

    1. Ensure that any molecules that have been broken across periodic
       bondaries have been made whole.
    2. Only supply molecules that :program:`propka` will support.

    To generate a timeseries of pKa predictions for a protein, first create
    a :class:`PropkaTraj` object by supplying an MDAnalysis Universe or
    AtomGroup::

      import MDAnalysis as mda
      # here we will use the MDAnalysis PSF and DCD example files
      from MDAnalysisTests.datafiles import PSF, DCD
      from propkatraj import PropkaTraj

      u = mda.Universe(PSF, DCD)

      pkatraj = PropkaTraj(u)

    By default :class:`PropkaTraj` will select all protein residues (based on
    standard protein residue naming) in the input ``atomgroup``, however this
    can be changed by passing the ``select`` keyword. For example to limit the
    input residues to residues numbered 1 through to 10::

      pkatraj = PropkaTraj(u, select="resnum 1-10")

    Alternatively, you can create an atomgroup beforehand and pass `None` to
    ``select``::

      ag = u.select_atoms('resnum 1-10")
      pkatraj = PropkaTraj(ag, select=None)

    Once the :class:`PropkaTraj` object has been created, we then use the
    :meth:`run` to analyse the trajectory::

      pkatraj.run()

    In some cases, especially considering that calling propka can be time
    consuming for large systems, you may not want to analyse every frame.
    You can pass the ``start``, ``stop``, and ``step`` arguments to :meth:`run`
    in order to alter the starting frame, final frame and frame sampling
    frequency. Please note that frame numbers are zero formatted, and the
    ``stop`` keyword is exclusive, as is the default in
    :class:`MDAnalysis.AnalysisBase` objects. For example, in order to sample
    every two frames from the first 10 frames::

      pkatraj.run(start=0, stop=10, step=2)

    You can also set the ``verbose`` keyword in order to have a nice progress
    bar of the analysis progress::

      pkatraj.run(verbose=True)

    Once completed, the resulting timeseries of predicted pKas per residue is
    made available as a :class:`pandas.DataFrame` under :attr`PropkaTraj`.pkas.
    For example if one wanted to get a summary of the statistics of the
    timeseries::

      pkatraj.pkas.describe()

    If one wanted to plot per residue boxplots of the timeseries data and save
    it as a file `pKa-plot.png` (note: this requires the ``matplotlib`` and
    ``seaborn`` packages which are not installed by default with
    ``propkatraj``)::

      import matplotlib.pyplot as plt
      import seaborn as sns

      fig = plt.figure(figsize=(28, 8))
      ax = fig.add_subplot(1, 1, 1)
      sns.boxplot(data=pkatraj.pkas, ax=ax)
      ax.set_xlabel('residue number')
      ax.set_ylabel(r'p$K_a$')
      fig.savefig('pKa-plot.png')


    .. versionadded:: 1.1.0
    """
    def __init__(self, atomgroup, select='protein', skip_failure=False,
                 **kwargs):
        if select is not None:
            self.ag = atomgroup.select_atoms(select).atoms
        else:
            self.ag = atomgroup.atoms

        # Issue #23 (keep until the PDBWriter is fixed)
        if len(self.ag.select_atoms('not protein')) > 0:
            wmsg = ("Non protein atoms passed to propka 3.1.\n MDAnalysis' "
                    "PDBWriter does not currently write non-standard residues "
                    "correctly as HETATM records and this may lead to "
                    "incorrect pKa predictions.\n"
                    "See https://github.com/Becksteinlab/propkatraj/issues/24 "
                    " for more details")
            warnings.warn(wmsg)

        self.tmpfile = os.path.join(os.path.curdir, 'current.pdb')
        self.skip_failure = skip_failure
        super(PropkaTraj, self).__init__(self.ag.universe.trajectory, **kwargs)

    def _prepare(self):
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
            # TODO: it would be nice to allow for other options, maybe for 3.2?
            mol = pk.single(pstream, optargs=['--quiet'])
        except (IndexError, AttributeError) as err:
            errmsg = "failure on frame: {0}".format(self._ts.frame)
            if not self.skip_failure:
                raise_from(RuntimeError(errmsg), err)
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
                self._columns = [g.atom.resNumb for g in groups]
        finally:
            # deallocate stream
            pstream.close(force=True)

    def _conclude(self):
        # Ouput failed frames
        if self.num_failed_frames > 0:
            perc_failure = (self.num_failed_frames / self.n_frames) * 100

            # if frames have failed we need to ammend times accordingly
            failed_indices = [np.where(self.times == i) for i in
                              self.failed_times]
            self.times = np.delete(self.times, failed_indices)

            wmsg = ("number of failed frames = {0}\n"
                    "percentage failure = {1}\n"
                    "failed frames: {2}".format(self.num_failed_frames,
                                                perc_failure,
                                                self.failed_frames_log))
            logging.warning(wmsg)

        self.pkas = pd.DataFrame(self._pkas, index=pd.Float64Index(self.times,
                                 name='time'), columns=self._columns)


@deprecate(release="1.1.0", remove="2.0.0",
           message="Use ``PropkaTraj(u, ..).run().pkas instead.")
def get_propka(universe, sel='protein', start=None, stop=None, step=None,
               skip_failure=False):
    """Get and store pKas for titrateable residues along trajectory.

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


    Notes
    -----
    Currently, temporary :program:`propka` files are written in the same
    directory as the input trajectory file. This will leave a ``current.pka``
    and ``current.propka_input`` file post-analysis. These are the temporary
    files for the final frame and can be removed. Should the trajectory file
    not have an input directory (e.g. when using MDAnalysis' `fetch_mmtf`
    method), then the files will be written to the current directory.

    Known issues:

    1. Due to the current behaviour of the MDAnalysis PDBWriter, non-protein
       atoms are written to PDBs using `ATOM` records instead of `HETATM`.
       This is likely to lead to undefined behaviour in :program:`propka`,
       which will likely expect `HETATM` inputs. We recommend users to only
       pass protein atoms for now. See the following issue for more details:
       https://github.com/Becksteinlab/propkatraj/issues/24

    """

    # need AtomGroup to write out for propka
    if isinstance(sel, string_types):
        atomsel = universe.select_atoms(sel)
    elif isinstance(sel, (list, np.ndarray)):
        atomsel = universe.atoms[sel]

    # Issue #23 (keep until the PDBWriter is fixed)
    if len(atomsel.select_atoms('not protein')) > 0:
        wmsg = ("Non protein atoms passed to propka 3.1.\n MDAnalysis' "
                "PDBWriter does not currently write non-standard residues "
                "correctly as HETATM records and this may lead to "
                "incorrect pKa predictions.\n"
                "See https://github.com/Becksteinlab/propkatraj/issues/24 "
                " for more details")
        warnings.warn(wmsg)

    # "filename" for our stream
    # use same name so that propka overwrites
    try:
        newname = os.path.join(os.path.dirname(universe.filename),
                               'current.pdb')
    except TypeError:
        # we have a trajectory without a directory
        newname = os.path.join(os.path.curdir, 'current.pdb')

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
            # https://github.com/Becksteinlab/propkatraj/issues/13
            # https://github.com/Becksteinlab/propkatraj/issues/10
            err_msg = "{0} (failure {2}): failing frame {1}".format(
                    universe.trajectory.filename, ts.frame, failed_frames)
            if not skip_failure:
                raise_from(RuntimeError(err_msg), err)
            else:
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
        logging.warning('percent failure = {0:.3f}%'.format(
            float(failed_frames)/len(universe.trajectory)*100))
        logging.warning('failed frames: %r', failed_frames_log)

    # a `pandas.DataFrame` is a good data structure for this data
    df = pd.DataFrame(pkas, index=pd.Float64Index(times, name='time'),
                      columns=[g.atom.resNumb for g in groups])

    return df
