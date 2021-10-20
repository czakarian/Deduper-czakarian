#!/usr/bin/env python

""" This program performs deduplicating of PCR duplicates from RNA-seq data. ... """

import argparse
import Bioinfo
from queue import PriorityQueue
import re 

def get_strand(flag:int) -> str:
    """This function returns the strand for a read in a SAM file based on the bit flag value. 
    Returns + if positive strand and - if reverse strand."""
    if(flag & 16) == 16:
        return "-"
    return "+"

def adjust_start_pos(start_pos:int, cigar:str) -> int:
    """This function adjusts for soft clipping and returns the corrected start position. If no soft clipping, will return back the inputted start_position."""
    pattern = '^([0-9]+)S'
    result = re.match(pattern, cigar)
    if result:  
        return start_pos - int(result.group(1))
    return start_pos

def get_position(line:str) -> str:
    """This function reads in a read line from a SAM file and generates a string with the position info (chrom, start position, strand)"""        
    cols = line.split()
    # extract the chromosome 
    chrom = cols[2] 
    if len(chrom) == 1:
        chrom = "0" + chrom
    # extract the start position (adjusted for soft clipping)  
    start_pos = str(adjust_start_pos(int(cols[3]), cols[5]))
    # modify the start position string to be of length 12 (add 0's to the beginning) // needed for correct sorting by PriorityQueue()
    start_pos = start_pos.rjust(12 - len(start_pos), '0')
    # extract the strand 
    strand = get_strand(int(cols[1]))
    out = chrom + ":" + str(start_pos) + ":" + strand
    return out

def get_args():
    """This function returns the parser arguments entered in command line"""
    parser = argparse.ArgumentParser(description="A program to remove PCR duplicates from SAM file")
    parser.add_argument("-f", "--file", help="Absolute file path of SAM file", required=True)
    parser.add_argument("-p", "--paired", help="Designates file is paired end", action='store_true')
    parser.add_argument("-u", "--umi", help="Designates file containing list of UMIs, unset if randomers used", nargs='?')
    return parser.parse_args()

# store the command line args in variables
args = get_args()
file= args.file
paired = args.paired
umi_file = args.umi

# Store list of known UMIs in a set 
known_UMIs = set()
with open(umi_file, "r") as fr:
    for line in fr:
        known_UMIs.add(line.strip())

# PriorityQueue object to store ...
# will use PriorityQueue() set to max length of 300 to store last 300 encountered positions [will be ordered by position] 
# every time we add a new position, will pop off the item in queue with lowest start position 
pos_queue = PriorityQueue(maxsize=300)

# Dictionary to store positions in queue at any specific moment 
# Key = Position, Value = set of UMIs
pos_UMI_dict = {}

with open(file, "r") as fr, open("test/test1.sam", "w") as fw:
    for line in fr:
        # if we are at header line, write line to output file
        if line[0] == "@":
            fw.write(line)
        else:
            cols = line.split()
            # extract and store the umi 
            umi = cols[0].split(":")[-1] 

            if umi in known_UMIs:
                pos = get_position(line)
                if pos in pos_UMI_dict:
                    if umi not in pos_UMI_dict[pos]:
                        pos_UMI_dict[pos].add(umi)
                        fw.write(line)
                else:
                    pos_UMI_dict[pos] = {umi}
                    if pos_queue.full():
                        pos_UMI_dict.pop(pos_queue.get())
                    pos_queue.put(pos)
                    fw.write(line)
