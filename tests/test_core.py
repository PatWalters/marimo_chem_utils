import pandas as pd
from marimo_chem_utils import smi2inchi_key, add_fingerprint_column

def test_smi2inchi_key():
    """
    Test the smi2inchi_key function.
    """
    assert smi2inchi_key("CCO") == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert smi2inchi_key("invalid_smiles") is None

def test_add_fingerprint_column():
    """
    Test the add_fingerprint_column function.
    """
    df = pd.DataFrame({"SMILES": ["CCO", "c1ccccc1"]})
    df = add_fingerprint_column(df, fp_type="fp")
    assert "fp" in df.columns
    assert df["fp"].iloc[0] is not None
    assert df["fp"].iloc[1] is not None
