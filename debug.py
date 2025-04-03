from Bio.PDB import PDBParser

pdb_file = "XP_500_bestposes_complex.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", pdb_file)

ligand_atoms = []
protein_atoms = []

for model in structure:
    for chain in model:
        for residue in chain:
            hetfield, resseq, icode = residue.id
            # Use standard ATOM records for protein (or residue id with blank hetfield)
            if residue.get_id()[0] == " " and residue.resname != "LIG":
                protein_atoms.extend(residue.get_atoms())
            elif residue.resname == "LIG":
                ligand_atoms.extend(residue.get_atoms())

print(f"Ligand atoms: {len(ligand_atoms)}")
print(f"Protein atoms: {len(protein_atoms)}")
