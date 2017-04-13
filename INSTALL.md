# Installation

## Pre-requisites

Install [MDAnalysis](http://mdanalysis.org),
[MDSynthesis](http://mdsynthesis.readthedocs.io),
[pandas](http://pandas.pydata.org/) and their dependencies:

    pip install pandas
    pip install mdanalysis mdanalysistests
    pip install mdsynthesis

[PROPKA 3.1](https://github.com/jensengroup/propka-3.1) does not have
official releases at the moment but we can
[pip-install](https://pip.pypa.io/en/stable/reference/pip_install/#vcs-support)
the source directly from GitHub (you need to have
[git](https://git-scm.com/) installed):

    pip install git+https://github.com/jensengroup/propka-3.1.git@master#egg=propka-3.1

This *should* install the most recent version of PROPKA 3.1 that
supports all required features for propkatraj.

## `propkatraj`

Copy `propkatraj.py` into your working directory or add the directory
containing this file to your Python `sys.path`.

The following should work:

    from propkatraj import get_propka

The `propkatra.get_propka()` function contains all functionality. 

