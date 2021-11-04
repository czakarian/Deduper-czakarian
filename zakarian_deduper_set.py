#!/usr/bin/env python

""" This program performs deduplicating/removal of PCR duplicates from single-end RNA-seq data. 
Requires as input a sorted SAM file of uniquely mapped reads and a list of UMIs. If list of UMIs
not provided, will assume random UMIs used."""

import argparse
# import Bioinfo 
import re 

def get_strand(flag:int) -> bool :
    """This function returns the strand for a read in a SAM file based on the bit flag value. 
    Returns True if positive strand and False if reverse strand."""
    return (flag & 16) != 16

def get_start(start_pos:int, cigar:str, strand:bool) -> int:
    """This function returns the true 5' start position (adjusted for soft clipping) 
    for reads mapping to + strand or the - strand."""
    
    # for + strand reads, adjust soft clipping at the beginning of CIGAR string
    if strand:
        result = re.match('^([0-9]+)S', cigar)
        if result:  
            start_pos -= int(result[1])
    
    # for - strand reads, sum M/D/N/S values and add that to start_pos, ignore any inserts and soft clipping at the beginning   
    else:
        result_S = re.search('([0-9]+)S$', cigar)
        result_MDN = re.findall('([0-9]+)[MDN]', cigar)
        if result_S:
            start_pos += int(result_S[1])
        if result_MDN:
            start_pos += sum(map(int, result_MDN))
        # subtract 1 to correct to actual position // technically not necessary for the purpose of identifying dups 
        start_pos -= 1
    
    return start_pos

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
invalid_UMIs = 0 # counter to track number of reads with invalid UMIs
dups = 0 # counter to track number of PCR duplicates removed
unique_reads = {"1":0} # counter dictionary for number of unique reads per chromosome
chrom = "1" # track the chromosome location while parsing
umi_set = set() # set of known umis (if provided)
umi_pos_strand_set = set() # Set of tuples of form (umi, start_pos, strand) to check for duplicates while parsing

# If UMI list provided, add umis to umi_set
if umi_file:
    with open(umi_file, "r") as fr:
        for line in fr:
            umi_set.add(line.strip())

with open(file, "r") as fr, open(file[:-4] + "_deduped.sam", "w") as fw:
    for line in fr:
        # write header lines to output file
        if line[0] == "@":
            fw.write(line)
        else:
            cols = line.split()

            # extract and store UMI, corrected start postion, and strand
            umi = cols[0].split(":")[-1] 
            strand = get_strand(int(cols[1]))
            start_pos = get_start(int(cols[3]), cols[5], strand)

            # if we reach a read on the next chromosome, store new chrom and clear set before parsing reads of next chrom 
            if chrom != cols[2]:
                chrom = cols[2]
                unique_reads[chrom] = 0
                umi_pos_strand_set = set()

            # Check if valid UMI: no N's and in umi_set if using known umi's
            if "N" not in umi and (not umi_file or (umi_file and umi in umi_set)):
                # if not a duplicate add to set and write to file 
                if (umi, start_pos, strand) not in umi_pos_strand_set:
                    umi_pos_strand_set.add((umi, start_pos, strand))
                    unique_reads[chrom] += 1 
                    fw.write(line)
                else:
                    dups += 1
            else:
                invalid_UMIs += 1


# Output summary statistics
print("Duplicates Removed: ", dups)
print("Invalid UMIs: ", invalid_UMIs)
print("Unique Reads Per Chromosome:")
for c in unique_reads:
    print(c, "\t", unique_reads[c])