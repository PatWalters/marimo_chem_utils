"""
A collection of cheminformatics utility functions for use with marimo notebooks.
"""
from typing import Optional, Tuple

import pandas as pd
import useful_rdkit_utils as uru
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
from tqdm.auto import tqdm
import altair as alt
import marimo as mo
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from rdkit.Chem import rdDepictor


def add_fingerprint_column(df: pd.DataFrame, column_name: str = None, fp_type: str = "counts_fp",
                           smiles_column: str = "SMILES") -> pd.DataFrame:
    """
    Adds a fingerprint column to a DataFrame.

    Args:
        df: DataFrame to add the fingerprint column to.
        column_name: Name of the new fingerprint column.
        fp_type: Type of fingerprint to generate.
                    Options: "fp", "counts_fp", "np_fp", "np_counts_fp".
        smiles_column: Name of the column containing SMILES strings.

    Returns:
        DataFrame with the added fingerprint column.
    """
    if column_name is None:
        column_name = fp_type
    tqdm.pandas()
    smi2fp = uru.Smi2Fp()
    fp_function_dict = {
        "fp": smi2fp.get_fp,
        "counts_fp": smi2fp.get_count_fp,
        "np_fp": smi2fp.get_fp,
        "np_counts_fp": smi2fp.get_np_counts
    }
    df[column_name] = df[smiles_column].progress_apply(fp_function_dict[fp_type])
    return df


def add_image_column(df: pd.DataFrame, image_column: str = "image", smiles_column: str = "SMILES") -> pd.DataFrame:
    """
    Adds a column with molecule images to a DataFrame.

    Args:
        df: DataFrame to add the image column to.
        image_column: Name of the new image column.
        smiles_column: Name of the column containing SMILES strings.

    Returns:
        DataFrame with the added image column.
    """
    tqdm.pandas()
    if image_column not in df.columns:
        image_list = df[smiles_column].progress_apply(uru.smi_to_base64_image, target='altair')
        df.insert(0, image_column, image_list)  # insert at the beginning of the dataframe
    return df


def smi2inchi_key(smiles: str) -> Optional[str]:
    """
    Converts a SMILES string to an InChIKey.

    Args:
        smiles: SMILES string to convert.

    Returns:
        InChIKey string, or None if the conversion fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToInchiKey(mol)
    return None


def draw_molecule_grid(df: pd.DataFrame,
                       smiles_column: str = "SMILES",
                       smarts=None,
                       template=None,
                       legend_column: str = None,
                       num_cols: int = 5,
                       image_size: Tuple[int, int] = (200, 200),
                       max_to_show: int = 25) -> Image.Image:
    """
    Draws a grid of molecules from a DataFrame.

    Args:
        df: DataFrame containing molecule data.
        smiles_column: Name of the column with SMILES strings.
        legend_column: Name of the column to use for legends under the images.
        num_cols: Number of columns in the grid.
        image_size: Size of each molecule image.
        max_to_show: Maximum number of molecules to display.

    Returns:
        A PIL Image object containing the grid.
    """
    rdDepictor.SetPreferCoordGen(True)
    df_to_show = df.head(max_to_show)
    smiles_list = df_to_show[smiles_column].tolist()
    # Handle structure alignment
    if template is not None:
        mol_list = uru.align_mols_to_template(template, smiles_list)
    else:
        mol_list = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]
    # Handle SMARTS highlighting
    match_list = []
    if smarts is not None:
        pat = Chem.MolFromSmarts(smarts)
        match_list = [mol.GetSubstructMatch(pat) if mol is not None else () for mol in mol_list]
    # Handle legends
    legend_list = []
    if legend_column is not None:
        if df_to_show[legend_column].dtype == float:
            legend_list = [f"{x:.2f}" for x in df_to_show[legend_column].tolist()]
        else:
            legend_list = [str(x) for x in df_to_show[legend_column]]
    mol_grid = Draw.MolsToGridImage(mol_list, molsPerRow=num_cols, subImgSize=image_size, legends=legend_list, highlightAtomLists=match_list)
    return mol_grid


def add_inchi_key_column(df: pd.DataFrame, inchi_key_column: str = "inchi_key",
                         smiles_column: str = "SMILES") -> pd.DataFrame:
    """
    Adds an InChIKey column to a DataFrame.

    Args:
        df: DataFrame to add the InChIKey column to.
        inchi_key_column: Name of the new InChIKey column.
        smiles_column: Name of the column containing SMILES strings.

    Returns:
        DataFrame with the added InChIKey column.
    """
    tqdm.pandas()
    df[inchi_key_column] = df[smiles_column].progress_apply(smi2inchi_key)
    return df


def add_tsne_columns(df: pd.DataFrame, smiles_column: str = "SMILES", fp_column: str = "np_count_fp") -> pd.DataFrame:
    """
    Adds TSNE_x and TSNE_y columns to a DataFrame.

    This function calculates fingerprints (if not already present),
    reduces their dimensionality using PCA, and then applies t-SNE
    to generate 2D coordinates.

    Args:
        df: DataFrame to add the t-SNE columns to.
        smiles_column: Name of the column containing SMILES strings.
        fp_column: Name of the column containing fingerprints.
                   If it doesn't exist, it will be calculated using 'np_counts_fp' type.

    Returns:
        DataFrame with the added TSNE_x and TSNE_y columns.
    """
    if fp_column not in df.columns:
        df = add_fingerprint_column(df, column_name=fp_column, fp_type="np_counts_fp", smiles_column=smiles_column)

    fp_array = np.vstack(df[fp_column].values)
    n_samples = fp_array.shape[0]

    if n_samples <= 1:
        df["TSNE_x"] = np.nan
        df["TSNE_y"] = np.nan
        return df

    # Reduce dimensions with PCA
    n_pca_components = min(50, n_samples - 1)
    pca = PCA(n_components=n_pca_components, random_state=42)
    fp_reduced = pca.fit_transform(fp_array)

    # Apply t-SNE
    # Perplexity must be less than n_samples
    perplexity = min(30.0, float(n_samples - 1))
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity, init='pca', learning_rate='auto')
    tsne_results = tsne.fit_transform(fp_reduced)

    df["TSNE_x"] = tsne_results[:, 0]
    df["TSNE_y"] = tsne_results[:, 1]

    return df


def interactive_chart(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: Optional[str] = None,
    cutoff_value: Optional[float] = None,
    x_title: str = "X",
    y_title: str = "Y",
    image_col: str = "image",
):
    """
    Creates an interactive Altair scatter plot with molecule tooltips and box selection.

    Args:
        df: The dataframe containing the data.
        x_col: The column name for the x-axis.
        y_col: The column name for the y-axis.
        color_col: The column to use for coloring points.
        cutoff_value: A cutoff to create a binary color scheme.
                      If None, a continuous color scale is used.
        x_title: The title for the x-axis.
        y_title: The title for the y-axis.
        image_col: Column with base64 images for tooltips.
    """
    # Define which columns to show in the tooltip
    tooltip_cols = [image_col]
    if color_col:
        tooltip_cols.append(color_col)

    # Base chart
    chart = alt.Chart(df).mark_circle(size=60)

    # Set up color encoding
    if color_col:
        if cutoff_value is not None:
            # Binary color scheme based on a cutoff
            color_encoding = alt.condition(
                alt.datum[color_col] > cutoff_value,
                alt.value("red"),
                alt.value("blue")
            )
        else:
            # Continuous color scale
            color_encoding = alt.Color(
                color_col,
                scale=alt.Scale(scheme='viridis'),
                legend=alt.Legend(title=color_col)
            )
    else:
        # Default color if no color column is provided
        color_encoding = alt.value("steelblue")

    # Build the final chart with encodings
    final_chart = chart.encode(
        x=alt.X(
            x_col,
            title=x_title,
            axis=alt.Axis(titleFontSize=16, labelFontSize=12)
        ),
        y=alt.Y(
            y_col,
            title=y_title,
            axis=alt.Axis(titleFontSize=16, labelFontSize=12)
        ),
        color=color_encoding,
        tooltip=[alt.Tooltip(c) for c in tooltip_cols]
    )

    return mo.ui.altair_chart(final_chart)

