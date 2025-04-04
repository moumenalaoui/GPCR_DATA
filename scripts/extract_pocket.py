import numpy as np
import copy
from Bio.PDB import PDBParser, PDBIO, Structure

# Settings
pdb_file = "XP_500_bestposes_complex.pdb"  # your merged complex file
output_file = "XP_500_bestposes_pocket.pdb"  # output file for the pocket
distance_threshold = 6.0  # in Å; adjust if needed

# Parse the structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", pdb_file)

# Collect ligand atoms (assume ligand residue name is "LIG")
ligand_atoms = []
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.get_resname() == "LIG":
                ligand_atoms.extend(list(residue.get_atoms()))

if not ligand_atoms:
    print("No ligand atoms found. Please check that your ligand residue is labeled 'LIG'.")
    exit()

# Identify protein residues that are within the threshold of any ligand atom
pocket_residues = []
for model in structure:
    for chain in model:
        for residue in chain:
            # Skip ligand residues
            if residue.get_resname() == "LIG":
                continue
            include_residue = False
            for atom in residue:
                for ligand_atom in ligand_atoms:
                    dist = np.linalg.norm(atom.coord - ligand_atom.coord)
                    if dist <= distance_threshold:
                        include_residue = True
                        break
                if include_residue:
                    break
            if include_residue:
                # We'll copy the residue to avoid modifying the original structure
                pocket_residues.append(copy.deepcopy(residue))

if not pocket_residues:
    print("No pocket residues found within the threshold.")
    exit()

# Create a new structure containing only the pocket residues.
# We'll build a new structure with one model and one chain.
from Bio.PDB import Model, Chain

new_structure = Structure.Structure("pocket")
new_model = Model.Model(0)
new_chain = Chain.Chain("A")
for residue in pocket_residues:
    new_chain.add(residue)
new_model.add(new_chain)
new_structure.add(new_model)

# Write the pocket structure to a PDB file
io = PDBIO()
io.set_structure(new_structure)
io.save(output_file)
print("Pocket saved as:", output_file)
