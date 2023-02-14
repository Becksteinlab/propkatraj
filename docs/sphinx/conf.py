# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

import propkatraj

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'propkatraj'
copyright = '2016-2023, David Dotson, Oliver Beckstein, Armin Zjajo, Rick Sexton, Shujie Fan, Irfan Alibay, Ian Kenney'
author = 'David Dotson'
release = '1.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autosummary',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Configuration for intersphinx: refer to the Python standard library
# and other packages used by MDAnalysis
intersphinx_mapping = {
                       'https://docs.python.org/3/': None,
                       'https://docs.mdanalysis.org/stable/index.html': None,
                       'https://numpy.org/doc/stable/': None,
                       'https://pandas.pydata.org/pandas-docs/stable/': None,
                       }
