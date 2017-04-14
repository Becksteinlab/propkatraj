# propkatrj --- https://github.com/Becksteinlab/propkatraj
# Copyright (c) 2013-2017 David Dotson, Armin Zjajo, Oliver Beckstein
# Released under the GNU General Public License v3+

from __future__ import print_function

from six import string_types
import os
import cStringIO

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis as mda


def get_propka(universe, sel='protein', start=None, stop=None, step=None):
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

    Results
    -------
    sim : :class:`mdsynthesis.Sim`
        Sim that was operated on. Returned here for convenience; useful with `dask`.

    """

    # need AtomGroup to write out for propka
    if isinstance(sel, string_types):
        atomsel = universe.select_atoms(sel)
    elif isinstance(sel, (list, np.array)):
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
    for ts in universe.trajectory[start:stop:step]:
        pm.echo(ts.frame, time=ts.time)

        # we create a named stream to write the atoms of interest into
        pstream = mda.lib.util.NamedStream(cStringIO.StringIO(), newname)
        atomsel.write(pstream)

        pstream.reset()         # reset for reading

        # we feed the stream to propka, and it reads it as if it were a file on
        # disk
        mol = pk.single(pstream, optargs=['--quiet'])
        pstream.close(force=True)  # deallocate

        # parse propka data structures to get out what we actually want
        confname = mol.conformation_names[0]
        conformation = mol.conformations[confname]
        groups = conformation.get_titratable_groups()

        # extract pka estimates from each residue
        pkas.append([g.pka_value for g in groups])

        # record time
        times.append(ts.time)

    # a `pandas.DataFrame` is a good data structure for this data
    df = pd.DataFrame(pkas, index=pd.Float64Index(times, name='time'),
                      columns=[g.atom.resNumb for g in groups])

    return df
