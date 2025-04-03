from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Read the ligand SMILES from the file
with open("ligand_fixed.smi", "r") as f:
    ligand_smiles = f.read().strip()
print("Ligand SMILES:", ligand_smiles)

# Try to create a molecule without sanitizing
mol = Chem.MolFromSmiles(ligand_smiles, sanitize=False)
if mol is None:
    raise ValueError("Mol could not be generated with sanitize=False.")

# Now attempt to sanitize manually, catching errors:
try:
    Chem.SanitizeMol(mol)
except Exception as e:
    print("Error during sanitization:", e)
    # If sanitization fails, you can try to fix aromaticity manually
    # For example, you might try:
    # Chem.SetAromaticity(mol)
    # But with such a large peptide, it may be necessary to use an SDF representation or MolStandardize tools.
    raise e

# Calculate descriptors
ligand_descriptors = {
    "mol_wt": Descriptors.MolWt(mol),
    "logp": Descriptors.MolLogP(mol),
    "num_h_donors": Descriptors.NumHDonors(mol),
    "num_h_acceptors": Descriptors.NumHAcceptors(mol)
}
ligand_df = pd.DataFrame([ligand_descriptors])

# Read the interaction features CSV
interaction_df = pd.read_csv("interaction_features.csv")

# Merge DataFrames
features_matrix = pd.concat([interaction_df, ligand_df], axis=1)

# Save the final feature matrix to CSV
features_matrix.to_csv("features_matrix.csv", index=False)
print("features_matrix.csv created!")
