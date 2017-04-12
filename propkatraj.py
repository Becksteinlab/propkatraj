from six import string_types
import cStringIO

import pandas as pd
import numpy as np

import propka.run as pk
import MDAnalysis.lib.util as util

def get_propka(sim, sel='protein', start=None, stop=None, step=1):
    """Get and store pKas for titrateable residues near the binding site.

    Parameters
    ----------
    sim : :class:`mdsynthesis.Sim`
        Sim to obtain pKas for.
    sel : str, array_like
        Selection string to use for selecting atoms to use from Sim's universe.
        Can also be a numpy array or list of atom indices to use.
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
        atomsel = sim.universe.select_atoms(sel)
    elif isinstance(sel, (list, np.array)):
        atomsel = sim.universe.atoms[sel]

    # "filename" for our stream
    newname = sim["current.pdb"].abspath  # use same name so that propka overwrites

    times = []
    pkas = []
    for ts in sim.universe.trajectory[start:stop:step]:
        print '\rTime (ps): {}'.format(ts.time), 
        
        # we create a named stream to write the atoms of interest into
        pstream = util.NamedStream(cStringIO.StringIO(), newname)
        atomsel.write(pstream)

        pstream.reset()         # reset for reading

        # we feed the stream to propka, and it reads it as if it were a file on
        # disk
        mol = pk.single(pstream, None)
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
    df = pd.DataFrame(pkas, index=pd.Float64Index(times, name='time (ps)'),
                      columns=[g.atom.resNumb for g in groups])

    # save DataFrame in HDF5 using `datreant.data`
    sim.data['propka/pka'] = df

    return sim
