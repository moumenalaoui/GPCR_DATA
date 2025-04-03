from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Read the ligand SMILES from the file (using the corrected file, e.g., ligand_fixed.smi)
with open("ligand_fixed.smi", "r") as f:
    ligand_smiles = f.read().strip()
print("Ligand SMILES:", ligand_smiles)

# Create an RDKit molecule without sanitization
mol = Chem.MolFromSmiles(ligand_smiles, sanitize=False)
if mol is None:
    raise ValueError("Mol could not be generated with sanitize=False.")

# Create sanitization flags that exclude kekulization
ops = Chem.SanitizeFlags.SANITIZE_ALL & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE
try:
    Chem.SanitizeMol(mol, sanitizeOps=ops)
except Exception as e:
    print("Error during sanitization (excluding kekulization):", e)
    raise e

# Optionally, if you need to clear aromatic flags (but for descriptor calculation it might be acceptable to leave them):
# Chem.Kekulize(mol, clearAromaticFlags=True)

# Calculate some basic descriptors for the ligand
ligand_descriptors = {
    "mol_wt": Descriptors.MolWt(mol),
    "logp": Descriptors.MolLogP(mol),
    "num_h_donors": Descriptors.NumHDonors(mol),
    "num_h_acceptors": Descriptors.NumHAcceptors(mol)
}
ligand_df = pd.DataFrame([ligand_descriptors])

# Read the interaction features CSV
interaction_df = pd.read_csv("interaction_features.csv")

# Merge the two DataFrames (assuming one sample per file)
features_matrix = pd.concat([interaction_df, ligand_df], axis=1)

# Save the final feature matrix to CSV
features_matrix.to_csv("features_matrix.csv", index=False)
print("features_matrix.csv created!")
