# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'JWSToolKit'
copyright = '2025, Valentin Delabrosse'
author = 'Valentin Delabrosse'
release = '1.0'

import os
import sys

sys.path.insert(0, os.path.abspath("../src"))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.autodoc", "sphinx.ext.viewcode", "sphinx.ext.napoleon"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo' #'sphinx_rtd_theme'
html_static_path = ['_static']


html_theme_options = {
    "light_logo" : 'JWSToolKit_logo_light_fullsize.png', 
    "dark_logo" : 'JWSToolKit_logo_dark_fullsize.png',
    "sidebar_hide_name": True  
}

