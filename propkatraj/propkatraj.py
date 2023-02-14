# propkatrj --- https://github.com/Becksteinlab/propkatraj
# Copyright (c) 2013-2020 David Dotson, Irfan Alibay, Ricky Sexton,
#                         Armin Zjajo, Shujie Fan, Oliver Beckstein
# Released under the GNU General Public License v2+

from io import StringIO
import os

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

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
    atomgroup : :class:`MDAnalysis.core.universe.Universe` or :class:`MDAnalysis.core.groups.AtomGroup`
        Group of atoms containing the residues for pKa analysis. Please note
        that :class:`MDAnalysis.core.groups.UpdatingAtomGroup` are *not* supported and will
        be automatically converted to :class:`MDAnalysis.core.groups.AtomGroup`.
    select : str
        Selection string to use for selecting a subsection of atoms to use
        from the input ``atomgroup``. Note: passing non-protein residues to
        :program:`propka` may lead to incorrect results (see notes).
        [`protein`]
    skip_failure : bool
        If set to ``True``, skip frames where :program:`propka` fails. A list
        of failed frames is made available in
        :attr:`PropkaTraj.results.failed_frames_log`. If ``False`` raise a
        :exc:`RuntimeError` exception on those frames. [`False`]
        
    Attributes
    ----------
    times : np.ndarray
        times of the successfully analyzed trajectory frames
    results.pkas : pd.DataFrame
        computed pKa's for each residue as a column and each frame as a row;
        the column names are the residue numbers
    results.num_failed_frames : int
        If PROPKA failed for any frames, contains number of failed frames.
        (Needs `skip_failure` to be ``True``.)
    results.failed_frames_log : list
        Frame indices of failed frames, if `skip_failure` set to ``True``
    results.failed_times : list
        Times of failed frames.


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
        super().__init__(self.ag.universe.trajectory, **kwargs)

    def _prepare(self):
        self._pkas = []
        self.results.pkas = []
        self.results.num_failed_frames = 0
        self.results.failed_frames_log = []
        self.results.failed_times = []
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
                raise RuntimeError(errmsg) from err
            else:
                warnings.warn(errmsg)
                self.results.num_failed_frames += 1
                self.results.failed_frames_log.append(self._ts.frame)
                self.results.failed_times.append(self._ts.time)
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
        if self.results.num_failed_frames > 0:
            perc_failure = (self.results.num_failed_frames / self.n_frames) * 100

            # if frames have failed we need to ammend times accordingly
            failed_indices = [np.where(self.times == i) for i in
                              self.results.failed_times]
            self.times = np.delete(self.times, failed_indices)

            wmsg = ("number of failed frames = {0}\n"
                    "percentage failure = {1}\n"
                    "failed frames: {2}".format(self.results.num_failed_frames,
                                                perc_failure,
                                                self.results.failed_frames_log))
            logging.warning(wmsg)

        self.results.pkas = pd.DataFrame(self._pkas, index=pd.Index(self.times,
                                         name='time', dtype='float64'),
                                         columns=self._columns)
