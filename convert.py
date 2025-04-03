from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# Read the SMILES from your file (assuming it has one SMILES string)
with open("ligand_fixed.smi", "r") as f:
    ligand_smiles = f.read().strip()

mol = Chem.MolFromSmiles(ligand_smiles)
ligand_descriptors = {
    "mol_wt": Descriptors.MolWt(mol),
    "logp": Descriptors.MolLogP(mol),
    "num_h_donors": Descriptors.NumHDonors(mol),
    "num_h_acceptors": Descriptors.NumHAcceptors(mol)
}
ligand_df = pd.DataFrame([ligand_descriptors])

# Read interaction features
interaction_df = pd.read_csv("interaction_features.csv")

# Combine (for one sample, you can simply concatenate columns)
features_matrix = pd.concat([interaction_df, ligand_df], axis=1)

# Save final feature matrix
features_matrix.to_csv("features_matrix.csv", index=False)
print("features_matrix.csv created!")
