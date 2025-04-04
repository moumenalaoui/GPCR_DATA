from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

# --- LIG_Feat: Compute descriptors for 5-HT using its SMILES ---
with open("5HT.smi", "r") as f:
    serotonin_smiles = f.read().strip()
print("5-HT SMILES:", serotonin_smiles)

# Create an RDKit molecule from the SMILES string
mol = Chem.MolFromSmiles(serotonin_smiles)
if mol is None:
    raise ValueError("Could not parse the 5-HT SMILES string.")

# Calculate standard small-molecule descriptors for 5-HT
ligand_descriptors = {
    "mol_wt": Descriptors.MolWt(mol),
    "logp": Descriptors.MolLogP(mol),
    "num_h_donors": Descriptors.NumHDonors(mol),
    "num_h_acceptors": Descriptors.NumHAcceptors(mol)
}
ligand_df = pd.DataFrame([ligand_descriptors])
print("Ligand descriptors:", ligand_df)

# --- INT_Feat: Load the interaction features CSV ---
interaction_df = pd.read_csv("interaction_features.csv")
print("Interaction features:", interaction_df)

# --- POCK_Feat: Compute a simple pocket descriptor (atom count) from the MOL2 file ---
def parse_mol2(file):
    with open(file, "r") as f:
        lines = f.readlines()
    # Find the start of the ATOM section
    atom_start = None
    for i, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            atom_start = i + 1
            break
    if atom_start is None:
        raise ValueError("No ATOM section found in the MOL2 file.")
    atom_lines = []
    for line in lines[atom_start:]:
        if line.startswith("@<TRIPOS>"):
            break
        if line.strip():
            atom_lines.append(line.strip())
    return atom_lines

atoms = parse_mol2("XP_500_bestposes_pocket.mol2")
pocket_atom_count = len(atoms)
pocket_descriptors = {"pocket_atom_count": pocket_atom_count}
pocket_df = pd.DataFrame([pocket_descriptors])
print("Pocket descriptors:", pocket_df)

# --- Combine All Features into One Feature Matrix ---
features_matrix = pd.concat([interaction_df, ligand_df, pocket_df], axis=1)
features_matrix.to_csv("features_matrix.csv", index=False)
print("features_matrix.csv created!")
