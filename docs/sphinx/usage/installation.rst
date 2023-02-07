Installation
============

Using ``pip``
-------------

The latest release version of propkatraj can be installed using *pip* and
requires a Python version of 3.8 or later.

.. code:: bash

    pip install propkatraj

From source
-----------

``propkatraj`` can also be installed from source:

.. code:: bash

    git clone https://github.com/Becksteinlab/propkatraj.git
    cd propkatraj
    pip install .

If installing this way, it is recommended to run the included tests.

Tests
~~~~~

After installing the base package, running tests require the additional
installation of the ``MDAnalysisTests`` and ``pytest`` packages. These can be
installed with either ``pip`` or ``conda``.

.. code:: bash

    // with pip
    pip install MDAnalysisTests pytest

    // with conda
    conda install MDAnalysisTests pytest

To run the tests, in the source root:

.. code:: bash

    pytest -v propkatraj/tests