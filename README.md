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
* [MDSynthesis](http://mdsynthesis.readthedocs.io)
* [pandas](http://pandas.pydata.org/)

## Documentation

See the Jupyter notebook
[docs/propkatraj-example.ipynb](./docs/propkatraj-example.ipynb) for
how to use `propkatraj.get_propka` on an example trajectory and how to
plot the data with [seaborn](https://seaborn.pydata.org/).

## Contact

Please raise issues in the
[issue tracker](https://github.com/Becksteinlab/propkatraj/issues).
