input_file = "XP_500_bestposes_5HT.pdb"
output_file = "XP_500_bestposes_5HT_repaired.pdb"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        # Only modify HETATM lines that are for the ligand (residue name "0")
        if line.startswith("HETATM"):
            # The residue name is in columns 18-20 (index 17:20)
            resname = line[17:20].strip()
            if resname == "0":
                # Replace with "LIG" (and keep the rest of the line as is)
                newline = line[:17] + "LIG" + line[20:]
                outfile.write(newline)
            else:
                outfile.write(line)
        else:
            outfile.write(line)

print("Repaired file saved as:", output_file)
