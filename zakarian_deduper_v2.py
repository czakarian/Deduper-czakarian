#!/usr/bin/env python

""" This program performs deduplicating/removal of PCR duplicates from single-end RNA-seq data. 
Requires as input a sorted SAM file of uniquely mapped reads and a list of UMIs. If list of UMIs
not provided, will assume random UMIs used."""

import argparse
import Bioinfo 

def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="This program removes PCR duplicates from a sorted SAM file of uniquely mapped reads.")
    parser.add_argument("-f", "--file", help="File path of SAM file", required=True)
    parser.add_argument("-p", "--paired", help="Designates a paired end SAM file", action='store_true')
    parser.add_argument("-u", "--umi", help="Designates file containing list of known UMIs, if not provided will assume random UMIs")
    return parser.parse_args()

# store the command line args in variables
args = get_args()
file= args.file
paired = args.paired
umi_file = args.umi

if paired:
    exit("Warning: This script does not yet handle paired end files. Exiting.")

# Initialize variables 
UMI_pos_dict = {} # Dictionary to store UMIs and encountered positions (Key = UMI, Value = Set of tuples of form: (start,strand) // dont need to store chromosome because dictionary will be reset for each chromosome
invalid_UMIs = 0 # counter to keep track of the total number of invalid UMIs encountered
dups = 0 # counter to keep track of how many PCR duplicates we encounter 
chrom = "1" # track which chromosome is being parsing so we don't have to store chrom individually for each read in our dictionary
unique_reads = {"1":0} # tracks # of unique reads per chromosome

# If UMI list provided, initialize as keys in dictionary
if umi_file:
    with open(umi_file, "r") as fr:
        for line in fr:
            UMI_pos_dict[line.strip()] = set() # set value as an empty set 

with open(file, "r") as fr, open(file[:-4] + "_deduped2.sam", "w") as fw:
    for line in fr:
        # write header lines to output file
        if line[0] == "@":
            fw.write(line)
        else:
            cols = line.split()
            # extract and store the umi 
            umi = cols[0].split(":")[-1] 

            # if we reach a read on the next chromosome, store new chrom value and clear out the dictionary before continuing 
            if cols[2] != chrom:
                chrom = cols[2]
                # if known UMIs, clear values of dictionary but keep UMIs as keys
                if umi_file:
                    UMI_pos_dict = {k: set() for k in UMI_pos_dict.keys()} # will reset values of all umi keys to empty set 
                # if random UMIs, just clear out the whole dictionary (keys and values)
                else:
                    UMI_pos_dict = {}
                unique_reads[chrom] = 0

            # extract and store the corrected start postion and strand 
            strand = Bioinfo.get_strand(int(cols[1]))
            start_pos = Bioinfo.get_start(int(cols[3]), cols[5], strand)

            # if umi in dict, its either (1) known UMI or (2) random UMI we have already encounterd 
            # for both cases, we will check if its a dup by checking if we have already seen a read at the position 
            if umi in UMI_pos_dict:
                # if (start_pos, strand) not in set for this umi, add to set and write read to output file 
                if (start_pos, strand) not in UMI_pos_dict[umi]:
                    UMI_pos_dict[umi].add((start_pos, strand))
                    fw.write(line)
                    unique_reads[chrom] += 1
                # otherwise it is a dup 
                else:
                    dups += 1

            # if umi not in dict and using known UMIs, it is an invalid UMI 
            # if using randomers, will mark as invalid if contains N
            elif umi_file or "N" in umi:
                invalid_UMIs += 1
            
            # if umi not in dict and using ramdomers, cant be dup so add umi to dictionary and write read to output
            else:
                UMI_pos_dict[umi] = {(start_pos, strand)}
                fw.write(line)
                unique_reads[chrom] += 1

# Output summary statistics
print("Duplicates Removed: ", dups)
print("Invalid UMIs: ", invalid_UMIs)
print("Unique Reads Per Chromosome:")
for c in unique_reads:
    print(c, "\t", unique_reads[c])