#! /usr/bin/python
"""Setuptools-based setup script for propkatraj.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages

setup(name='propkatraj',
      version='0.1.0-dev',
      description='obtain pKas for titreatable residues from a simulation trajectory',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        ],
      packages=['propkatraj.py'],
      scripts=[],
      license='GPLv3',
      long_description=open('README.md').read(),
      install_requires=['numpy', 'pandas', 'MDAnalysis', 'propka']
      )
