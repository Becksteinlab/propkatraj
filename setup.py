#! /usr/bin/python
"""Setuptools-based setup script for propkatraj.

For a basic installation just type the command::

  python setup.py install

"""

import sys
from setuptools import setup, find_packages
import versioneer

# Make sure we have the right python version (2.7+)
if sys.version_info[:2] < (2, 7):
    print('propkatraj requires python 2.7 or above')
    sys.exit(-1)

setup(name='propkatraj',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='obtain pKas for titreatable residues from a simulation trajectory',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      classifiers=['Development Status :: 6 - Mature',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',
                   'Programming Language :: Python :: 3.8',
      ],
      project_urls={
          'Documentation': 'https://github.com/Becksteinlab/propkatraj/blob/master/README.md',
          'Source': 'https://github.com/Becksteinlab/propkatraj',
          'Issue Tracker': 'https://github.com/Becksteinlab/propkatraj/issues',
      },
      packages=['propkatraj'],
      scripts=[],
      license='GPLv3',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown; variant=GFM',
      install_requires=['six', 'numpy', 'pandas', 'MDAnalysis<2.0.0', 'propka==3.1']
      )
