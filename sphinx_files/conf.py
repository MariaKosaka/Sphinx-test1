# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, '/Users/maria/Desktop/Labo/B4M1Seminar-2023/tutorial-20230710/Sphinx_test/src/')

project = 'Sphinx Study'
copyright = '2023, maria'
author = 'maria'
release = '1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
 'sphinx.ext.autodoc',
 'sphinx.ext.viewcode',
 'sphinx.ext.todo',
 'sphinx.ext.napoleon'
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'ja'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]


#LOGO
html_logo = './_static/LOGO1.PNG'
#html_theme_options = {
#'logo_only': True, 
#'prev_next_buttons_location': 'bottom'
#}
