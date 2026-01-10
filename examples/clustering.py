# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.19.1",
#     "marimo-chem-utils==0.2.2",
#     "pandas==2.3.3",
#     "rdkit==2025.9.3",
#     "tqdm==4.67.1",
#     "useful-rdkit-utils==0.93",
# ]
# ///

import marimo

__generated_with = "0.19.1"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Building a Simple Cluster Viewer with Marimo

    Clustering provides an efficient means of understand chemcial datasets.  However, while a number of clustering algorithms such as Butina and BitBirch are available, there aren't many tools available that enable us to view the clustering results.  In this notebook we'll look at how we can use [marimo](https://marimo.io/) to quickly build a viewer that will enable us to interactively view clustering results.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 1. Import the Necessary Python Libraries
    We begin by importing a few key Python libraries
    """)
    return


@app.cell
def _():
    import pandas as pd
    import marimo_chem_utils as mcu
    import useful_rdkit_utils as uru
    from rdkit import Chem
    from rdkit import rdBase
    from tqdm.auto import tqdm
    import marimo as mo
    return Chem, mcu, mo, pd, rdBase, tqdm, uru


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 2. Read the Input Data, Add Fingerprints, and Cluster
    Next, we read the input data.  In this case the data is a set of carbonic anyhydrase inhibitors taken from a 2018 paper by [Isidro Cortés-Ciriano and Andreas Bender](https://pubs.acs.org/doi/10.1021/acs.jcim.8b00542).  After reading the data, we use the [marimo_chem_utils](https://github.com/PatWalters/marimo_chem_utils) library to add a fingerprint column to the dataframe, then use the [useful_rdkit_utils](https://github.com/PatWalters/useful_rdkit_utils) library to cluster the chemical structures using the Butina alogorithm as implemneted in the RDKit.
    """)
    return


@app.cell
def _(mcu, pd, uru):
    df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/datafiles/refs/heads/main/carbonic.csv")
    _ = mcu.add_fingerprint_column(df,fp_type="fp")
    df['cluster'] = uru.taylor_butina_clustering(df.fp)
    return (df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 3. Align the Structures for Each Cluster

    When examining closely related chemical structures, it's always easier when the structures are aligned.  In this step, we loop over clusers then use the function [mcs_align](https://useful-rdkit-utils.readthedocs.io/en/latest/misc_utils.html#useful_rdkit_utils.misc_utils.mcs_align) from useful_rdkit_utils to align the structures.  Note this procedure uses the RDKit's [AllChem.GenerateDepictionMatching2DStructure](https://www.herongyang.com/Cheminformatics/RDKit-GenerateDepictionMatching2DStructure-Orientation.html) to generate the aligned structures.  This method works well, but if the aligned structrue is very long along the horizontal or vertical dimension, the depiction algorithm may rotate it a bit.
    """)
    return


@app.cell
def _(Chem, df, pd, rdBase, tqdm, uru):
    cluster_df_list = []
    cluster_rep_list = []
    for k,v in tqdm(df.groupby("cluster")):    
        if len(v) > 1:
            with rdBase.BlockLogs():
                v['mol'] = uru.mcs_align(v.SMILES)
        else:
            v['mol'] = v.SMILES.apply(Chem.MolFromSmiles)
        cluster_df_list.append(v)
        cluster_rep_image = uru.mol_to_base64_image(v.mol.values[0],target="altair")
        cluster_rep_list.append([k,cluster_rep_image, len(v)])
    aligned_df = pd.concat(cluster_df_list)
    cluster_rep_df = pd.DataFrame(cluster_rep_list,columns=["cluster","mol","count"])
    return aligned_df, cluster_rep_df


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 4. Set Up the "Controller" Table
    The viewer we construct below will have two panes.  The pane on the left has one representative from each cluster. Clicking on the check box in the far left column, selects a cluster and the first (up to 12) members of that cluster are shown in the pane on the right.  To set up this "controller" table, we will use marimo's [mo.ui.table](https://docs.marimo.io/api/inputs/table/) class.  We set a few parameters to improve the workflow.
    - selection="single" - Only allows the table to select single rows
    - page_size=5 - set to ensure that the structures are large enough to see.  Note that tooltips automatically provide a larger version of the structure.
    - show_column_summaries=False - Turns off the summary graphs in the table header
    - show_data_types=False - Turns of the display of datatypes in the header
    """)
    return


@app.cell
def _(cluster_rep_df, mo):
    cluster_rep_table = mo.ui.table(cluster_rep_df,selection="single",page_size=5,show_column_summaries=False,show_data_types=False)
    return (cluster_rep_table,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 5. Set Up the "Detail Table"

    Next we set up  the chemical structure view that is shown on the right in the view below.  If nothing is selected in the controller table, we use marimo's markdown component [mo.md]() to put up a message asking the user to make a selection in the table on the right. If a row is selected, we collect the RDKit molecules corresponding to the selected cluster and use the RDKit's [Chem.Draw.MolsToGridImage](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html#rdkit.Chem.Draw.MolsToGridImage) to display them. Note that we can get the selected cluster by accessing `cluster_rep_table.value` which holds a Pandas dataframe containing the selected row.
    """)
    return


@app.cell
def _(Chem, aligned_df, cluster_rep_table, mo):
    if len(cluster_rep_table.value):
        selected_cluster_id = cluster_rep_table.value.cluster.values[0]
        cluster_grid_df = aligned_df.query(f"cluster == {selected_cluster_id}").head(12)
        legends = cluster_grid_df.Name.astype("str").tolist()
        cluster_grid = Chem.Draw.MolsToGridImage(cluster_grid_df.mol.values.tolist(),molsPerRow=3,legends=legends)
    else:
        cluster_grid = mo.md("Please make a selection from the table on the left")
    return (cluster_grid,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### 6. Combine the Two Views

    Now that we've generated the controller and detail views on our data, we can simply stack them in a horizontal view using marimo's [mo.hstack](https://docs.marimo.io/api/layouts/stacks/#marimo.hstack).  To use the hhview, we simply pass it a list containing the two views and second list containing the preferred proportionla widths. Note that you can click on the "Fullscreen" icon to the right of the view to zoom it to full screen.
    """)
    return


@app.cell(hide_code=True)
def _(cluster_grid, cluster_rep_table, mo):
    mo.hstack([cluster_rep_table,cluster_grid],widths=[1,5])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
