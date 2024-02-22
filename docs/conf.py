# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath("."))

# -- Project information -----------------------------------------------------

project = "Yet Another Viral Subspecies Analysis Pipeline"
copyright = "2022, Thomas A. Christensen II"
author = "Thomas A. Christensen II"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "paramparse",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.githubpages",
    "sphinx.ext.intersphinx",
    "sphinx_multiversion",
]

myst_enable_extensions = [
    "strikethrough",
    "deflist",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "_ext", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_logo = "images/yavsap_logo.png"
html_theme_options = {
    "style_external_links": True,
}

html_context = {
    "display_github": True,
    "github_repo": "ksumngs/yavsap",
}

github_url = "https://github.com/ksumngs/v-met/yavsap"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ["_static"]

html_sidebars = {
    "**": [
        "versions.html",
    ],
}

smv_branch_whitelist = r".*(develop|master)"

intersphinx_mapping = {
    "nextflow": ("https://nextflow.io/docs/latest/", None),
    "singularity": ("https://apptainer.org/user-docs/master", None),
}
