# Installation

## Pre-requisites

Install [MDAnalysis](http://mdanalysis.org),
[pandas](http://pandas.pydata.org/) and their dependencies:

    pip install pandas
    pip install mdanalysis mdanalysistests

[PROPKA 3.1](https://github.com/jensengroup/propka-3.1) does not have
official releases at the moment but we can
[pip-install](https://pip.pypa.io/en/stable/reference/pip_install/#vcs-support)
the source directly from GitHub (you need to have
[git](https://git-scm.com/) installed):

    pip install git+https://github.com/jensengroup/propka-3.1.git@master#egg=propka-3.1

This *should* install the most recent version of PROPKA 3.1 that
supports all required features for propkatraj.

## `propkatraj`

Install from source with:

    git clone https://github.com/Becksteinlab/propkatraj.git
    cd propkatraj
    pip install .

Use the `--user` flag to install among your local user packages.

The following should then work:

    from propkatraj import get_propka

The `propkatra.get_propka()` function contains all functionality. 

