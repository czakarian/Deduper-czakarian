#!/usr/bin/env python

""" This program performs deduplicating of PCR duplicates from RNA-seq data. For the script to work, 
need to provide a sorted SAM file as input. Random UMIs or known UMIs. Single or Paired end. """

import argparse
#import Bioinfo
import re 

def get_strand(flag:int) -> bool :
    """This function returns the strand for a read in a SAM file based on the bit flag value. 
    Returns True if positive strand and False if reverse strand."""
    return (flag & 16) != 16

def get_start(start_pos:int, cigar:str, strand:bool) -> int:
    """This function returns the corrected start position (adjusted for soft clipping) and the actual start position for reads mapping to reverse strand."""
    adjust = 0
    if strand:
        result = re.match('^([0-9]+)S', cigar)
        if result:  
            adjust = 0 - int(result[1])
    else:
        # for - strand, sum M/D/N/S values and add that to start_pos, ignore inserts and soft clipping at the beginning   
        result_S = re.search('([0-9]+)S$', cigar)
        result_M = re.findall('([0-9]+)M', cigar)
        result_D = re.findall('([0-9]+)D', cigar)
        result_N = re.findall('([0-9]+)N', cigar)
        if result_S:
            adjust += int(result_S[1])
        if result_M:
            adjust += sum(map(int, result_M))
        if result_D:
            adjust += sum(map(int, result_D))
        if result_N:
            adjust += sum(map(int, result_N))
        # subtract 1 to correct for position 
        adjust -= 1
    return start_pos + adjust

def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="A program to remove PCR duplicates from SAM file")
    parser.add_argument("-f", "--file", help="Absolute file path of SAM file", required=True)
    parser.add_argument("-p", "--paired", help="Designates file is paired end", action='store_true')
    parser.add_argument("-u", "--umi", help="Designates file containing list of UMIs, unset if randomers used")
    return parser.parse_args()

# store the command line args in variables
args = get_args()
file= args.file
paired = args.paired
umi_file = args.umi

if args.paired:
    exit("Warning: This script does not yet handle paired end files. Exiting.\n")

# Initialize variables 
UMI_pos_dict = {} # Dictionary to store encountered positions under each UMI (Key = UMI, Value = Set of tuples storing start position and strand) // will be reset for each chromosome
invalid_UMIs = 0 # counter to keep track of the total number of invalid UMIs encountered
dups = 0 # counter to keep track of how many PCR duplicates we encounter 
chrom = "1" # will keep track of which chromosome we are currently on while parsing so we don't have to store chrom individually for each read 

count = 0
y = 0

# Store list of given UMIs in a set (if given) and initialize as keys in our dictionary
if umi_file:
    with open(umi_file, "r") as fr:
        for line in fr:
            line = line.strip()
            UMI_pos_dict[line] = set() # set value as an empty set
  

# changes for random UMIs?
# will add to dict if no N's, so add to the else part 


with open(file, "r") as fr, open("out.sam", "w") as fw:
    for line in fr:
        # if we are at header line, write line to output file
        if line[0] == "@":
            fw.write(line)
        else:
            cols = line.split()
            # extract and store the umi 
            umi = cols[0].split(":")[-1] 
            # if umi in given umi list proceed, otherwise skip read
            if umi in UMI_pos_dict:
                # if we reach next chromosome, clear out the dictionary before proceeding
                if cols[2] != chrom:
                    print(chrom, " ", count)
                    chrom = cols[2]
                    count = 0
                    UMI_pos_dict = {k: set() for k in UMI_pos_dict.keys()} # will reset values of all umi keys to empty set 
                # extract and store the starting postion and strand 
                strand = get_strand(int(cols[1]))
                start_pos = get_start(int(cols[3]), cols[5], strand)
                # if we haven't already encountered this pos/strand with this umi, add to set for that umi and write read to output file 
                if (start_pos, strand) not in UMI_pos_dict[umi]:
                    UMI_pos_dict[umi].add((start_pos, strand))
                    fw.write(line)
                    count += 1
                else:
                    if chrom == "Y":
                        y += 1
                    dups += 1
            else:
                invalid_UMIs += 1

print(chrom, " ", count)
print("Duplicates Removed: ", dups)
print("Invalid UMIs: ", invalid_UMIs)
print("Y", y)