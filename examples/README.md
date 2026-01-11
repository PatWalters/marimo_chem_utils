# marimo-Chem-utils Examples

This directory contains example [marimo](https://marimo.io/) notebooks demonstrating the usage of the `marimo_chem_utils` library for cheminformatics tasks.

## Examples

### 1. Clustering Viewer (`clustering.py`)
This notebook demonstrates how to build an interactive cluster viewer.
-   **Features:**
    -   Reads a dataset of carbonic anhydrase inhibitors.
    -   Generates fingerprints and performs Taylor-Butina clustering.
    -   Aligns structures within clusters using MCS alignment.
    -   Provides an interactive interface to select clusters and view aligned chemical structures.

### 2. REOS Filter Viewer (`reos.py`)
This notebook shows how to run and visualize Rapid Elimination of Swill (REOS) filters.
-   **Features:**
    -   Applies structural filters (e.g., BMS rules) to remove undesirable compounds.
    -   Summarizes filter results.
    -   Interactive table to select specific filter rules and view matching molecules with the matching substructure highlighted.

### 3. Interactive Scatterplot (`scatterplot.py`)
This notebook creates an interactive scatterplot for chemical space visualization.
-   **Features:**
    -   Calculates TSNE coordinates from chemical fingerprints.
    -   Displays an interactive scatterplot where points represent molecules.
    -   Allows users to select points on the plot to view the corresponding chemical structures in a grid.

## Usage

To run any of these examples, ensure you have the required dependencies installed (listed in the script metadata) and run the script using `marimo`:

```bash
marimo edit clustering.py
# or
marimo edit reos.py
# or
marimo edit scatterplot.py
```
