# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SCRIMP'
copyright = '2015-2022 Ghislain Haine'
author = 'Ghislain Haine'
release = '0.5.0'
version = '0.5.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import sys, os
sys.path.append(os.path.dirname(__file__))
sys.path.append('..')
sys.path.append(os.path.join('..','..'))
sys.path.append(os.path.join('..','scrimp'))
sys.path.append(os.path.join('..','scrimp','utils'))
sys.path.append(os.path.join('..','scrimp','examples'))

import sphinx_rtd_theme

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

autodoc_mock_imports = [
    'gmsh',
    'getfem',
    'matplotlib',
    'numpy',
    'scipy',
    'petsc4py',
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_favicon = 'Logo-ico.png'
html_logo = 'Logo.png'
html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'collapse_navigation': False,
    'navigation_depth': 4,
}

