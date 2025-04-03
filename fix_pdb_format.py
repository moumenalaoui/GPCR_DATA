input_file = "XP_500_bestposes_5HT_fixed.pdb"
output_file = "XP_500_bestposes_5HT_cleaned.pdb"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if line.startswith("HETATM"):
            # Fix residue name (17–20), chain (21), resSeq (22–26)
            new_line = (
                line[:17] + "LIG" + line[20] + "A" + "{:>4}".format("999") + line[26:]
            )
            outfile.write(new_line)
        else:
            outfile.write(line)

print("Fixed file saved as:", output_file)
