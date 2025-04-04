input_file = "XP_500_bestposes_5HT_cleaned.pdb"
output_file = "XP_500_bestposes_5HT_fixed_protein.pdb"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        # Check if line starts with HETATM and does NOT belong to the ligand (i.e. not LIG)
        if line.startswith("HETATM") and "LIG" not in line[17:20]:
            # Replace the record name from "HETATM" to "ATOM  "
            newline = "ATOM  " + line[6:]
            outfile.write(newline)
        else:
            outfile.write(line)

print("Fixed file saved as:", output_file)
