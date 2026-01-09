# /// script
# requires-python = ">=3.8"
# dependencies = [
#   "marimo",
#   "pandas==2.3.3",
#   "rdkit==2025.9.3",
#   "useful-rdkit-utils==0.96",
#   "marimo-chem-utils==0.2.4"
# ]
# ///

import marimo

__generated_with = "0.19.1"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## REOS in a [marimo](https://marimo.io) Notebook
    This notebook provides a quick overview of how the [Rapid Elimination of Swill (REOS)](https://practicalcheminformatics.blogspot.com/2018/08/filtering-chemical-libraries.html) filters can be run in a Marimo notebook.  We'll take a look at how the [useful_rdkit_utils](https://patwalters.github.io/Useful-RDKit-Utils/) library and marimo's [table](https://docs.marimo.io/api/inputs/table/) capability can be used to create a quick interactive viewer.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 1. Import the Necessary Python Libraries
    """)
    return


@app.cell
def _():
    from rdkit import Chem
    import useful_rdkit_utils as uru
    import marimo_chem_utils as mcu
    import pandas as pd
    import marimo as mo
    return Chem, mcu, mo, pd, uru


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2. Define a Support Fuctions
    Define a quick function to return the SMILES for the largest fragment in a molecule.
    """)
    return


@app.cell
def _(Chem, uru):
    def strip_salts(smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol =  uru.get_largest_fragment(mol)
            return Chem.MolToSmiles(mol)
        else:
            return None
    return (strip_salts,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3. Read and process a SMILES file.
    1. Read the SMILES file into a Pandas dataframe.
    2. Use the function we defined above to remove counterions, waters of hydration, etc.
    """)
    return


@app.cell
def _(pd, strip_salts):
    url = "https://raw.githubusercontent.com/PatWalters/datafiles/refs/heads/main/REOS/reos_test_1k.smi"
    df = pd.read_csv(url,sep=" ",names=["SMILES","Name"])
    df.SMILES = df.SMILES.apply(strip_salts)
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 4. Instantiate a REOS object.

    We will use the BMS rules with our REOS object.
    """)
    return


@app.cell
def _(uru):
    reos = uru.REOS(["BMS"])
    return (reos,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 5. Filter the Structures
    Run the filters on the dataframe `df`.  This applies the REOS SMARTS patterns to each of the SMILES in the dataframe.
    """)
    return


@app.cell
def _(df, reos):
    reos_df = reos.pandas_smiles(df.SMILES)
    df['reos'] = reos_df.description
    df['rule_set_name'] = reos_df.rule_set_name
    return (reos_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 6. Create a Summary Dataframe
    Create a new dataframe with a summary of the REOS run.
    """)
    return


@app.cell
def _(reos, reos_df):
    reos_summary_df = reos_df.value_counts(["rule_set_name","description"]).to_frame().reset_index()
    reos_dict = reos.get_rules_dict()
    reos_summary_df['SMARTS'] = [reos_dict.get((x,y)) for x,y in reos_summary_df[["rule_set_name","description"]].values]
    return (reos_summary_df,)


@app.cell
def _(reos_summary_df):
    reos_summary_df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 7. Set Up the Interactive Viewer
    With `reos_summary_df` in hand, we can use Marimo's `table` component to create an interactive viewer.  First, we'll pass the list of REOS filters that were triggered by the dataset and use this to create a table we can use to make selections.  To limit selections to only one row at a time, we set the `selection` parameter for the table to `single`.
    """)
    return


@app.cell
def _(mo, reos_summary_df):
    cols = ["description","count"]
    left_table = mo.ui.table(reos_summary_df[cols],selection="single",show_column_summaries=False,show_data_types=False)
    return (left_table,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The next cell contains the logic that coordinates the selection in the table on the left with the structure images shown on the right.  Once a selection is made in the table, the variable `left_table.value` changes and triggers the generation of a new chemical structure grid, which is shown to the right of the table.
    """)
    return


@app.cell
def _(df, left_table, mcu, mo, reos_summary_df):
    if len(left_table.value):
        selected_rule = left_table.value.description.values[0]
        grid_df = df.query(f"reos == '{selected_rule}'").head(12)
        reos_smarts = reos_summary_df.query(f"description == '{selected_rule}'").SMARTS.values[0]
        right_grid = mcu.draw_molecule_grid(grid_df,num_cols=4,smarts=reos_smarts)
    else:
        right_grid = mo.md("Please make a selection from the table on the left")
    return (right_grid,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 8. View the Results
    This is an interactive grid for viewing the results of the REOS calculation.  Clicking on a checkbox in the table on the left shows the structures that matched the query in the panel on the right.  The substructure matching the REOS rule is highlighted.
    """)
    return


@app.cell
def _(left_table, mo, right_grid):
    mo.hstack([left_table,right_grid],widths=[1,4])
    return


if __name__ == "__main__":
    app.run()
