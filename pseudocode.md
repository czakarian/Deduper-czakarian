# Deduper

### Problem + objective:
We want to develop a strategy to perform reference-based PCR duplicate removal for RNA-seq data. PCR duplicates are identical reads derived from sequencing of fragments that originated from the same original unamplified fragment. We do not want to include these reads in our transcript counts since they do not represent a real transcript. In order to differentiate PCR duplicates from transcripts that coincidentally fragmented the same way, we will use unique molecular identifiers (UMIs). Since each transcript should have a unique UMI, any identical reads that end up having the same UMI can be identified as being PCR duplicates. Reads that are identical (have the same start position) but differ in UMIs can each be counted as real transcripts. Our code will remove duplicates from single-end read data taking into account UMIs and soft-clipping. 

### Input:   
(1) a sorted SAM file with uniquely mapped reads 
(2) a text file with a list of 96 known UMIs
### Output:  
a SAM file with PCR duplicate reads removed

### SAM file columns of interest:
1. chromosome // col 3 // RNAME
2. position // col 4 // POS
3. strand // col 2 // FLAG // bit 16 set if reverse strand
4. soft clipping // col 6 // CIGAR 
5. UMI // col 1 // QNAME // split by ':' -> last element in list = UMI
 
### Pseudocode:
```
samtools sort # by first sorting by position, we wont have to check for pcr dups across all reads

Include an argparse option to specify text file with list of UMIs [if no arg given, will assume random UMIs used]
If using known UMIs, read in the UMIs from the provided text file and store in a set (known_UMIs)

Open the sorted input SAM file for reading and the output SAM file for writing 
    ->ouput @ header lines to the output SAM file

# my general strategy for parsing the file will be to use a PriorityQueue() object to store the last 300 encountered positions alongside a dictionary of the same size that will store the positions in queue +  encountered UMIs at each position
# this way we only need to store pos/UMI info for up to 300 items at a time instead of keeping an ongoing list of all pos/UMIs  
# with the queue every time we add a new element, we will also remove an element [with lowest position], allowing us to maintain at most 300 elements in the queue and dictionary at any one moment
# i am assuming that after adjusting for soft-clipping the adjusted position will never be >300 nucleotides away (impossible if our read lengths are less than 300)
# the max size of 300 can be adjusted if necessary, depending on if we know the max read length of the file, can set it to that value 

before starting to loop through reads, set up containers to store positions and UMIs
    pos_queue -> will use PriorityQueue() set to max length of 300 to store last 300 encountered positions [will be ordered by position] 
        every time we add a new position, will pop off the item in queue with lowest start position 
    pos_UMI_dict -> will store the positions in queue at any specific moment with the position as the key and the value a set of UMIs

for each line:
    read in the line & split by tabs [store in var called cols]
    extract the UMI from the QNAME (cols[0]) and store in var called cur_UMI
    if cur_UMI in known_UMIs: 
        extract & store the position using get_position() function which will take as input the read line and return a string of form 11:123456:+ [chr:start:strand] 
            store in var called cur_pos
        if cur_pos in pos_UMI_dict:
            if cur_UMI not in pos_UMI_dict[cur_pos]:
                add cur_UMI as element in the pos_UMI_dict[cur_pos] set
                output the read   
        else:
            add cur_pos to pos_UMI_dict with cur_UMI as an element in the set 
            if pos_queue is full (has 300 elements)
                pop off the element in pos_queue with lowest position + delete it from the dictionary
            add cur_pos to pos_queue 
            output the read 

Once we reach the end of the SAM file, we will close the input and output SAM files
```
### High Level Functions:

```
def get_strand(seq:int) -> str:
"""This function returns the strand(+/-) for a read in a SAM file based on the bit flag value"""
    # returns + if positive strand and - if reverse complemented strand 
    # takes in the bit flag value as input and checks whether the 16th bit is set (using & operator), - if set, + otherwise
    return a string of either + or - 
Example:
    Input: 36 
    Expected output: + 
    Input: 16 
    Expected output: - 
```
```
def adjust_start_pos(start_pos:str, cigar:str) -> int:
"""This function returns the corrected start position for a read that has been softclipped"""
    # Takes as input the start position and the CIGAR string to check for soft clipping 
    # Use regex to check for soft clipping at the beginning of cigar string ([0-9]+)S.
        if a match is found (# followed by an S), we want to substract the # from the start_pos
        else just return the same start_pos 
    return an integer representing the corrected start position 
Example:
    Input: 76814284, 2S71M 
    Expected output:  76814282
    Input: 76814284, 71M 
    Expected output:  76814284
```
```
def get_position(line:int) -> str:
"""This function reads in a read line from a SAM file and generates a string with the position info (chrom, start position, strand)"""        
    # split line by tab delimeter
    # grab the chromosome (3rd item in list)
    # call adjust_start_pos() on the start position (4th item in list)
    # call get_strand() on the bit flag (2nd item in list) to get + or - 
    return a string of the form ch:123456:s where ch is the chromosome, 123456 is the start position, and s is the strand (+ or -)
Example:
    Input: 
    Expected output: 11:123456:+
```