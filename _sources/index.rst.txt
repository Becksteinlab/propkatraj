.. propkatraj documentation master file, created by
   sphinx-quickstart on Mon Feb  6 10:22:22 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

propkatraj: ensamble pKa estimates
==================================

:Release: |release|
:Date: |today|
:Citation: |zenodo|

The :mod:`propkatraj` module contains an analysis class :class:`PropkaTraj`
which uses an ensemble approach for computing pKa estimates. It uses the single-structure
pKa estimator `propka`_ v3.1.0 across frames of a trajectory using the `MDAnalysis`_
package. `propkatraj`_ is hosted on GitHub and is released under the GPLv2
license or later version.

.. _propka: https://github.com/jensengroup/propka
.. _MDAnalysis: https://www.mdanalysis.org/
.. _propkatraj: https://github.com/Becksteinlab/propkatraj

.. |zenodo| image:: https://zenodo.org/badge/88095629.svg
   :alt: Zenodo DOI
   :target: https://zenodo.org/badge/latestdoi/88095629

.. toctree::
   :maxdepth: 4
   :hidden:

   usage/installation
   usage/examples
   api
