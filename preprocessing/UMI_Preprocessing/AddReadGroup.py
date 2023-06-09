#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# We load the BAM file and add the UG tag, a combination of the CB and UB tag.
# Then we remove all UMIs per cell that have less than --min_reads_per_umi reads.
# The number of reads per UMI was calculated in ReadMoleculeInfo.R.

# Libraries.
import pysam
import sys
import pandas as pd

# Loading the supplied input.
# sys.argv[0] is the script itself.
# sys.argv[1] the input bam file.
# sys.argv[2] the output folder.
# sys.argv[3] the UG file.
bam_input_path = sys.argv[1]
output_path = sys.argv[2]
UGs_path = sys.argv[3]
print("Using:")
print("BAM:    " + bam_input_path)
print("Output: " + output_path)
print("UGs:    " + UGs_path)

# We read the UG file.
ugs = pd.read_csv(UGs_path, header = None, names = ["UG"])

# We open the original BAM file.
bam_file = pysam.AlignmentFile(bam_input_path, "rb")

# We open a connection to the output file.
out_file = pysam.AlignmentFile(output_path+"possorted_genome_UG_tagged.bam", "wb", template=bam_file)

# We read the BAM file.
for line in bam_file:
    # line = next(bam_file) # Reading the next line. This is for test purposes only.
    # print(line)
    # We check if the CB tag is present.
    CB_check = line.has_tag("CB")
    # We check if the UB tag is present.
    UB_check = line.has_tag("UB")
    # We only proceed if both tags are present. Otherwise, we go to the next read.
    if all((CB_check, UB_check)):
        # We add the new tag.
        line.tags += [("UG", line.get_tag("CB")+"_"+line.get_tag("UB"))]
        # We check if the tag is white listed.
        if line.get_tag("UG") in ugs["UG"].values:
            # Writing the new line to the output file.
            out_file.write(line)

bam_file.close()
out_file.close()
