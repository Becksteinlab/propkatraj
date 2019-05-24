# README: propkatraj
[![DOI](https://zenodo.org/badge/88095629.svg)](https://zenodo.org/badge/latestdoi/88095629)

`propkatraj.py` can be used to computationally estimate pKa values for
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
* [MDAnalysis](https://mdanalysis.org)
* [pandas](https://pandas.pydata.org/)

See
[INSTALL.md](https://github.com/Becksteinlab/propkatraj/blob/master/INSTALL.md)
for how to install everything.

## Usage

The `propkatraj.get_propka()` function contains all
functionality. Import it with

```python
from propkatraj import get_propka
```

It takes a `MDAnalysis.Universe` instance as an argument and runs PROPKA on each
frame of the trajectory.

```
get_propka(universe, sel='protein', start=None, stop=None, step=None)

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
```

The function returns a
[pandas.DataFrame](http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe)
which contains the time as the first column and the residue numbers as
subsequent columns. For each time step, the predicted pKa value for
this residue is stored. Process the `DataFrame` to obtain statistics
as shown in the [Documentation](#Documentation).


## Documentation

See the Jupyter notebook
[docs/propkatraj-example.ipynb](https://nbviewer.jupyter.org/github/Becksteinlab/propkatraj/blob/master/docs/propkatraj-example.ipynb)
for how to use `propkatraj.get_propka` on an example trajectory and
how to plot the data with [seaborn](https://seaborn.pydata.org/).

## Citation

If you use `propkatraj` in published work please cite Reference 1 for
PROPKA 3.1 and Reference 2 for the ensemble method itself. Reference 3
is for the software if you need a specific software citation.

1. C. R. Søndergaard, M. H. M. Olsson, M. Rostkowski, and
   J. H. Jensen. Improved treatment of ligands and coupling effects in
   empirical calculation and rationalization of pKa values. *J
   Chemical Theory and Computation*, 7(7):2284–2295, 2011. doi:
   [10.1021/ct200133y](https://doi.org/10.1021/ct200133y).
   
2. C. Lee, S. Yashiro, D. L. Dotson, P. Uzdavinys, S. Iwata,
   M. S. P. Sansom, C. von Ballmoos, O. Beckstein, D. Drew, and
   A. D. Cameron. Crystal structure of the sodium-proton antiporter
   NhaA dimer and new mechanistic insights. *J Gen Physiol*,
   144(6):529–544, 2014. doi:
   [10.1085/jgp.201411219](https://doi.org/10.1085/jgp.201411219).

3. Oliver Beckstein, David Dotson, Rick Sexton, Shujie Fan, and Armin Zijajo. 
   (2019, May 24). Becksteinlab/propkatraj: 1.0.0 (Version release-1.0.0). 
   Zenodo. http://doi.org/10.5281/zenodo.3228426

## Contact

Please raise issues in the
[issue tracker](https://github.com/Becksteinlab/propkatraj/issues).
