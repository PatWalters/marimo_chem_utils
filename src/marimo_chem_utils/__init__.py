"""
A collection of cheminformatics utility functions for use with marimo notebooks.
"""
from importlib.metadata import version
from .core import (
    add_fingerprint_column,
    add_image_column,
    smi2inchi_key,
    draw_molecule_grid,
    add_inchi_key_column,
    add_tsne_columns,
    interactive_chart,
)

__version__ = version("marimo-chem-utils")

__all__ = [
    "add_fingerprint_column",
    "add_image_column",
    "smi2inchi_key",
    "draw_molecule_grid",
    "add_inchi_key_column",
    "add_tsne_columns",
    "interactive_chart",
]