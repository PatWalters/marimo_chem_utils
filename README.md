# marimo-chem-utils

A collection of cheminformatics utility functions for use with [marimo](https://marimo.io/), an open-source reactive notebook for Python.

## Installation

You can install `marimo-chem-utils` using `pip`:

```bash
pip install marimo-chem-utils
```

## Usage

Here is an example of how to use some of the functions in this library.

```python
import pandas as pd
import marimo as mo
from marimo_chem_utils import (
    add_fingerprint_column,
    add_image_column,
    add_inchi_key_column,
    add_tsne_columns,
    interactive_chart
)

# Load some data
df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/datafiles/refs/heads/main/carbonic.csv")

# Add fingerprints, images, and InChI keys
df = add_fingerprint_column(df, fp_type="counts_fp")
df = add_image_column(df, smiles_column="SMILES")
df = add_inchi_key_column(df, smiles_column="SMILES")

# Generate t-SNE coordinates
df = add_tsne_columns(df)

# Create an interactive plot
chart = interactive_chart(df, x_col="TSNE_x", y_col="TSNE_y", color_col="PIC50", image_col="image")

# In a marimo cell, you would display the chart like this:
# mo.ui.altair_chart(chart)
```

This will create a `pandas.DataFrame` with additional columns for fingerprints, images, and InChI keys, and then generate an interactive t-SNE plot that you can explore in your marimo notebook. The plot will display molecule images in the tooltips.

## Examples

The `examples` directory contains marimo notebooks that demonstrate how to use `marimo-chem-utils`. You can run them from the command line.

- **scatterplot.py**: An interactive scatterplot that shows chemical structures as tooltips. Selections on the plot are used to display a grid of molecules.

```bash
marimo edit examples/scatterplot.py
```

- **clustering.py**: An example of clustering molecules and viewing the cluster members.

```bash
marimo edit examples/clustering.py
```

