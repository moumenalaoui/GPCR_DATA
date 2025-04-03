from Bio.PDB import PDBParser
from scipy.spatial import distance
import pandas as pd
import os

pdb_file = "XP_500_bestposes_complex.pdb"
# Setup
parser = PDBParser(QUIET=True)
try:
    structure = parser.get_structure("complex", pdb_file)
except Exception as e:
    print(f"Error parsing PDB file: {e}")
    exit()
# Separate ligand and protein atoms
ligand_atoms = []
protein_atoms = []

for model in structure:
    for chain in model:
        for residue in chain:
            hetfield, resseq, icode = residue.id
            if hetfield.startswith("H") or residue.resname == "LIG":
                for atom in residue:
                    ligand_atoms.append(atom)
            elif hetfield == " ":
                for atom in residue:
                    protein_atoms.append(atom)

# Compute interaction types
close_contacts = 0
h_bond_like = 0
hydrophobic_contacts = 0

for l_atom in ligand_atoms:
    l_coord = l_atom.coord
    for p_atom in protein_atoms:
        p_coord = p_atom.coord
        dist = distance.euclidean(l_coord, p_coord)
        if dist <= 4.0:
            close_contacts += 1
            # H-bond like: O/N - O/N atoms within 3.5Å
            if (l_atom.element in ['O', 'N']) and (p_atom.element in ['O', 'N']) and dist <= 3.5:
                h_bond_like += 1
            # Hydrophobic proxy: C-C within 4Å
            if (l_atom.element == 'C') and (p_atom.element == 'C'):
                hydrophobic_contacts += 1

# Save result
df = pd.DataFrame([{
    "file_name": os.path.basename(pdb_file),
    "close_contacts": close_contacts,
    "h_bond_like": h_bond_like,
    "hydrophobic_contacts": hydrophobic_contacts
}])
df.to_csv("interaction_features.csv", index=False)
print("Interaction features saved to interaction_features.csv")
