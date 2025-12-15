import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Displaying Interactive Scatterplots with Chemical Structures
    In this notebook, we show how `marimo` can be used to create an interactive scatterplot that shows chemical structures as tooltips.  In addition, the plot allows selections which can be used to show the corresponding structures as a grid.  This capability is based on functions in the `marimo_chem_utils` library..
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Import the necessary Python libraries
    """)
    return


@app.cell
def _():
    import pandas as pd
    import marimo_chem_utils as mcu
    import marimo as mo
    return mcu, mo, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Read data from a csv file and sort by the pIC50 column.  The data is from a 2018 paper by [Isidro Cortés-Ciriano and Andreas Bender](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00542).
    """)
    return


@app.cell
def _(pd):
    df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/datafiles/refs/heads/main/carbonic.csv")
    df.sort_values("pIC50",ascending=False,inplace=True)
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Use the function `add_image_column` from `marimo_chem_utils` to add a chemical structure image to the dataframe. Note that a larger version of the structure appears as a tooltip when you hold the mouse over the image in the table.
    """)
    return


@app.cell
def _(df, mcu):
    if "image" not in df.columns:
        mcu.add_image_column(df)
    return


@app.cell
def _(df):
    df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Use the function `add_tnse_columns` from `marimo_chem_utils` to add two columns 'TSNE_x' and 'TSNE_y' to the datframe.  These columns are caculated by:
    - Generating a Morgan count fingerprint using the RDKit
    - Using PCA to reduce the fingerprint to 50 dimensions
    - Using TSNE to reduce the 50-dimensional vectors from PCA to 2 dimensions
    """)
    return


@app.cell
def _(df, mcu):
    _ = mcu.add_tsne_columns(df)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Generate an interactive scatterplot
    """)
    return


@app.cell
def _(df, mcu):
    tsne_chart = mcu.interactive_chart(df,"TSNE_x","TSNE_y",color_col="pIC50")
    return (tsne_chart,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Get the selected structures from the scatterplot and covert them to a molecule grid.
    """)
    return


@app.cell
def _(mcu, mo, tsne_chart):
    selected = tsne_chart.value
    if len(selected):
        mol_image = mcu.draw_molecule_grid(selected,legend_column="pIC50")
    else:
        mol_image = mo.md("Please make a selection")
    return (mol_image,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Display the molecule grid
    """)
    return


@app.cell
def _(mo, mol_image, tsne_chart):
    mo.hstack([tsne_chart,mol_image],widths="equal")
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
