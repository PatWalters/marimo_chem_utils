import marimo

__generated_with = "0.18.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import pandas as pd
    import useful_rdkit_utils as uru
    import marimo_chem_utils as mcu
    import marimo as mo
    return mcu, mo, pd, uru


@app.cell
def _(pd):
    df = pd.read_csv("https://raw.githubusercontent.com/PatWalters/datafiles/refs/heads/main/carbonic.csv")
    return (df,)


@app.cell
def _(df, mcu):
    mcu.add_fingerprint_column(df,fp_type="fp")
    return


@app.cell
def _(df, uru):
    df['cluster'] = uru.taylor_butina_clustering(df.fp)
    return


@app.cell
def _(df):
    df
    return


@app.cell
def _(df):
    cluster_rep_df = df.sort_values('cluster').drop_duplicates("cluster")[["SMILES","cluster"]]
    cluster_rep_df
    return (cluster_rep_df,)


@app.cell
def _(df, uru):
    cluster_count_df = uru.value_counts_df(df,"cluster")
    cluster_count_df
    return (cluster_count_df,)


@app.cell
def _(cluster_count_df, cluster_rep_df):
    cluster_df = cluster_rep_df.merge(cluster_count_df,on="cluster")
    cluster_df
    return (cluster_df,)


@app.cell
def _(cluster_df):
    cluster_df
    return


@app.cell
def _(cluster_df, mcu):
    mcu.add_image_column(cluster_df)
    cluster_df
    return


@app.cell
def _(cluster_df, mo):
    table = mo.ui.table(cluster_df[["image","cluster","count"]],selection="single",show_column_summaries=False)
    return (table,)


@app.cell
def _(df, mcu, mo, table):
    if len(table.value):
        cluster_id = table.value.cluster.values[0]
        selected_df = df.query(f"cluster == {cluster_id}").sort_values("pIC50",ascending=False)
        mol_image = mcu.draw_molecule_grid(selected_df,num_cols=2,legend_column="pIC50",max_to_show=10)
    else:
        mol_image = mo.md("Please click on a checkbox in the table")
    return (mol_image,)


@app.cell
def _(mo, mol_image, table):
    mo.hstack([table,mol_image],widths=[1,2])
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
