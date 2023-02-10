# Installation

## Pre-requisites

Install [MDAnalysis](https://mdanalysis.org),
[pandas](http://pandas.pydata.org/) and their dependencies, using
your preferred method such as using `conda`

    conda install mdanalysis mdanalysistests pandas
	
or 	`pip`:

    pip install pandas
    pip install mdanalysis mdanalysistests


## `pip` installation

Install `propkatraj` from [PyPi:
propkatraj](https://pypi.org/project/propkatraj/) with

    pip install propkatraj
	
which will install all additional dependencies such as [PROPKA
3](https://github.com/jensengroup/propka).


## Source installation

Install from source with:

    git clone https://github.com/Becksteinlab/propkatraj.git
    cd propkatraj
    pip install .

Use the `--user` flag for `pip` to install among your local user packages.

# Testing

Regression tests are provided under `propkatraj/tests`.

To run these you will need access to the
[pytest](https://docs.pytest.org/en/latest/index.html) package. This can be
installed either using `pip` or `conda`.

You can run the tests in the following way:

    pytest -v --disable-pytest-warnings propkatraj/tests

