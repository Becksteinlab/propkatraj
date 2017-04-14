# README: propkatraj

`propkatraj.py` can be used computationally estimate pKa values for
protein residues. We use an ensemble approach where many different
conformations are sampled with equilibrium molecular dynamics
simulations. We then apply the fast heuristic pKa predictor
[PROPKA 3.1](https://github.com/jensengroup/propka-3.1) to individual
frames of the trajectory. By analysing the statistics of the pKa
predictions a more consistent picture emerges than from a pKa
prediction of a single static conformation.


## Required software

* [PROPKA 3.1](https://github.com/jensengroup/propka-3.1) (used as a
  Python package)
* [MDAnalysis](http://mdanalysis.org)
* [pandas](http://pandas.pydata.org/)

See [INSTALL.md](INSTALL.md) for how to install everything.

## Usage

The `propkatra.get_propka()` function contains all functionality. 

    from propkatraj import get_propka

It takes a MDAnalysis.Universe instance as argument and runs PROPKA on each
frame of the trajectory.


    def get_propka(universe, sel='protein', start=None, stop=None, step=None):
        Get and store pKas for titrateable residues near the binding site.
    
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
        pkas : :class:`pandas.DataFrame`
            DataFrame giving estimated pKa value for each residue for each
            trajectory frame. Residue numbers are given as column labels, times as
            row labels.


## Documentation

See the Jupyter notebook
[docs/propkatraj-example.ipynb](./docs/propkatraj-example.ipynb) for
how to use `propkatraj.get_propka` on an example trajectory and how to
plot the data with [seaborn](https://seaborn.pydata.org/).

## Contact

Please raise issues in the
[issue tracker](https://github.com/Becksteinlab/propkatraj/issues).
