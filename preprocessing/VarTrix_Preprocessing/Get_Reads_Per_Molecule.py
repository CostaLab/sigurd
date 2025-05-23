#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

We load the BAM file and add the molecule tag.
It consists of the cell barcode (CB) and the UMI (UB).
We also count the number of occurences per molecule to further subset the BAM file.

"""

# We load the BAM file and add the UG tag. We also get the number of reads per molecule.

# Libraries.
import pysam
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Loading the supplied input.
bam_input_path = sys.argv[1]
output_path = sys.argv[2]
output_prefix = sys.argv[3]
sample = sys.argv[4]

print("Using:")
print("BAM:           " + bam_input_path)
print("Output:        " + output_path)
print("Output Prefix: " + output_prefix)
print("Sample:        " + sample)

# We open the original BAM file.
bam_file = pysam.AlignmentFile(bam_input_path, "rb")

# We open a connection to the output file.
out_file = pysam.AlignmentFile(output_path+output_prefix+".bam", "wb", template = bam_file)

# Dictionary to count the occurences of UG.
ug_counts = defaultdict(int)

# We read the BAM file.
for line in bam_file:
    # lines_read = lines_read + 1
    # print("Lines: "+str(lines_read))
    # line = next(bam_file) # Reading the next line.
    # print(line)
    # We only proceed if both tags are present. Otherwise, we go to the next read.
    if line.has_tag("CB") and line.has_tag("UB"):
        cb_tag = line.get_tag("CB")
        ub_tag = line.get_tag("UB")
        new_molecule = cb_tag + "_" + ub_tag
        line.set_tag("UG", new_molecule, value_type = "Z")
        out_file.write(line)
        ug_counts[new_molecule] += 1

bam_file.close()
out_file.close()

# We save the number of reads per molecule.
ug_counts_df = pd.DataFrame(list(ug_counts.items()), columns = ["Molecule", "Count"])
ug_counts_df[["CB", "UMI"]] = ug_counts_df["Molecule"].str.split("_", expand = True)
ug_counts_df["Rank"] = ug_counts_df["Count"].rank(ascending = False, method = "first")
ug_counts_df = ug_counts_df[["Molecule", "CB", "UMI", "Rank", "Count"]]
ug_counts_df = ug_counts_df.sort_values("Rank")

plt.figure(figsize=(5, 4))
plt.scatter(ug_counts_df["Rank"], ug_counts_df["Count"])
plt.xlabel("Rank")
plt.ylabel("Count")
plt.title("Molecule Rank vs. Count\n"+sample)
plt.savefig(output_path+output_prefix+"_MoleculeRank_"+sample+".png")

# We get the number of cells retained per molecule threshold.
thresholds = list(range(1, 11)) + list(range(20, 101, 10))
x_positions = np.arange(len(thresholds))
cells = []

for threshold in thresholds:
    count = ug_counts_df[ug_counts_df["Count"] >= threshold]["CB"].nunique()
    cells.append(count)

cells = pd.DataFrame({"Threshold": thresholds, "Cells": cells})
ug_counts_df["Rank"] = ug_counts_df["Count"].rank(ascending = False, method = "first")

plt.figure(figsize=(5, 4))
plt.plot(x_positions, cells["Cells"], marker="o")
plt.xticks(x_positions, thresholds)
plt.xlabel("Threshold")
plt.ylabel("Cells")
plt.title("Molecule Threshold vs. Cells\n"+sample)
plt.savefig(output_path+output_prefix+"_MoleculeCells_"+sample+".png")

ug_counts_df.to_csv(output_path+output_prefix+".csv", index = False)
