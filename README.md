# README: propkatraj
[![DOI](https://zenodo.org/badge/88095629.svg)](https://zenodo.org/badge/latestdoi/88095629)
[![GH Actions CI](https://github.com/Becksteinlab/propkatraj/actions/workflows/gh-ci.yaml/badge.svg?branch=main)](https://github.com/Becksteinlab/propkatraj/actions/workflows/gh-ci.yaml)
[![codecov](https://codecov.io/gh/Becksteinlab/propkatraj/branch/main/graph/badge.svg)](https://codecov.io/gh/Becksteinlab/propkatraj/branch/main)
[![docs](https://github.com/Becksteinlab/propkatraj/actions/workflows/docs.yaml/badge.svg?branch=main)](https://becksteinlab.github.io/propkatraj/)
[![Powered by MDAnalysis](https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA)](https://www.mdanalysis.org)

`propkatraj.py` can be used to computationally estimate pKa values for
protein residues. We use an ensemble approach where many different
conformations are sampled with equilibrium molecular dynamics
simulations. We then apply the fast heuristic pKa predictor
[PROPKA 3](https://github.com/jensengroup/propka) to individual
frames of the trajectory. By analysing the statistics of the pKa
predictions a more consistent picture emerges than from a pKa
prediction of a single static conformation.


## Required software

* [PROPKA 3](https://github.com/jensengroup/propka) (used as a
  Python package)
* [MDAnalysis](https://mdanalysis.org)
* [pandas](https://pandas.pydata.org/)

See
[INSTALL.md](https://github.com/Becksteinlab/propkatraj/blob/main/INSTALL.md)
for how to install everything.

## Usage

The `propkatraj.PropkaTraj` class contains all
functionality. Import it with

```python
from propkatraj import PropkaTraj
```

It takes a `MDAnalysis.AtomGroup` or `MDAnalysis.Universe` instance as an
argument to initialize and runs PROPKA on each frame of the trajectory when
calling the `run()` method. See `help(PropkaTraj)` for more details.

```python
pkatraj = PropkaTraj(atomgroup, select='protein', skip_failure=False)

   Runs :program:`propka` on the titrateable residues of the selected AtomGroup
   on each frame in the trajectory.
   
   Parameters
   ----------
   atomgroup : :class:`MDAnalysis.Universe` or :class:`MDAnalysis.AtomGroup`
       Group of atoms containing the residues for pKa analysis. Please note
       that :class:`MDAnalysis.UpdatingAtomGroup` are not supported and will
       be automatically converted to :class:`MDAnalysis.AtomGroup`.
   select : str
       Selection string to use for selecting a subsection of atoms to use
       from the input ``atomgroup``. Note: passing non-protein residues to
       :program:`propka` may lead to incorrect results (see notes). [`protein`]
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


pkatraj.run()

   Perform the calculation

   Parameters
   ----------
   start : int, optional
      start frame of analysis
   stop : int, optional
      stop frame of analysis
   step : int, optional
      number of frames to skip between each analysed frame
   verbose : bool, optional
      Turn on verbosity

```

Calling the `run()` method creates a [pandas.DataFrame](http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe),
accessed through `results.pkas`, which contains the time as the first column
and the residue numbers as subsequent columns. For each time step, the
predicted pKa value for this residue is stored. Process the `DataFrame` to
obtain statistics as shown in the [Documentation](#Documentation). For example,
you can get a summary of the statistics of the timeseries in the following
manner:

```python
pkatraj.results.pkas.describe()
```

## Documentation

See the Jupyter notebook
[docs/propkatraj-example.ipynb](https://nbviewer.jupyter.org/github/Becksteinlab/propkatraj/blob/main/docs/propkatraj-example.ipynb)
for how to use `propkatraj.PropkaTraj` on an example trajectory and
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

3. David Dotson, Irfan Alibay, Rick Sexton, Shujie Fan, Armin Zijajo, Oliver Beckstein. 
   (2020). Becksteinlab/propkatraj: 1.1.x. Zenodo. https://doi.org/10.5281/zenodo.3228425 

## Contact

Please raise issues in the
[issue tracker](https://github.com/Becksteinlab/propkatraj/issues).
