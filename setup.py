#! /usr/bin/python
"""Setuptools-based setup script for propkatraj.

For a basic installation just type the command::

  python setup.py install

"""

from setuptools import setup, find_packages
import versioneer

setup(name='propkatraj',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='obtain pKas for titreatable residues from a simulation trajectory',
      author='David Dotson',
      author_email='dotsdl@gmail.com',
      classifiers=['Development Status :: 6 - Mature',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   'Programming Language :: Python :: 3.8',
                   'Programming Language :: Python :: 3.9',
                   'Programming Language :: Python :: 3.10',
                   'Programming Language :: Python :: 3.11',
                   ],
      project_urls={
          'Documentation': 'https://github.com/Becksteinlab/propkatraj/blob/main/README.md',
          'Source': 'https://github.com/Becksteinlab/propkatraj',
          'Issue Tracker': 'https://github.com/Becksteinlab/propkatraj/issues',
          },
      packages=find_packages(),
      scripts=[],
      license='GPLv2+',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown; variant=GFM',
      python_requires='>=3.8',
      install_requires=['numpy', 'pandas', 'MDAnalysis>=2.0.0', 'propka==3.1']
      )
