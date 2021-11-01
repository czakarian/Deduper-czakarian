# Bioinfo module 
# Last update 11-1-21

import re

DNAbases = set('ATGCatcg')
RNAbases = set('AUGCaucg')
base_pairs = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}

def validate_base_seq(seq: str,RNAflag: bool =False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

def convert_phred(letter: str) -> int:
    """Converts a single character into a phred score"""
    return ord(letter)-33

def qual_score(phred_score: str) -> float:
    """Calculates the average quality score given a string of quality scores"""
    qsum = 0
    for letter in phred_score:
        qsum += convert_phred(letter)
    return qsum/len(phred_score)

def gc_content(DNA: str) -> float:
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA), "String contains invalid characters"
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

def calc_N50(ls: list) -> int:
    """This function takes as input a list of contig lengths and returns the N50 value of the assembly"""
    midpoint = 0.5 * sum(ls) # store the midpoint of the assembly
    ls_sorted = sorted(ls,reverse=True) # sort lengths in order from largest to smallest
    s = 0 # will store the sum of lengths until s >= midpoint 
    for i in ls_sorted:  
        s += i 
        if s >= midpoint:
            return i
    return 0

def oneLineFasta(input: str, output: str) -> None:
    """This function takes a FASTA file with multiple sequence lines per read and
    returns a FASTA file with one sequence line per read."""
    with open(input, "r") as fr, open(output, "w") as fw:
        first_header = True
        for line in fr:
            # if we are at a header line, output the header to file
            if line[0:1] == ">": 
                # make sure to add a new line char before each header (expect the first) since we strip the new line from end of the prev seq line
                if first_header != True:
                    line = "\n" + line
                first_header = False         
                fw.write(line)
            # if we are at a seq line, write to the file but strip new line char
            else:
                fw.write(line.strip())

def rev_comp(seq:str) -> str:
    """This function returns the reverse complement of a DNA sequence"""
    # Base pair dictionary to allow generation of a complement sequence 
    # Key = base, Value = complement base

    # generate the complement sequence
    comp_seq = ""
    for c in seq:
        comp_seq += base_pairs[c]
    # return reverse complemented sequence
    return comp_seq[::-1]

def meets_Qcutoff(qvalues: str) -> bool:
    """This function takes as input a string of index quality scores
    and returns True if none of the indexes have a quality score that
    is less than 30.""" 
    for c in qvalues:
        if convert_phred(c) < 30:
            return False
    return True

def get_strand(flag:int) -> bool :
    """This function returns the strand for a read in a SAM file based on the bit flag value. 
    Returns True if positive strand and False if reverse strand."""
    return (flag & 16) != 16

def get_start(start_pos:int, cigar:str, strand:str) -> int:
    """This function returns the corrected start position (adjusted for soft clipping) 
    and the actual 5' start position for reads mapping to reverse strand."""
    adjust = 0 # will track adjustment to make to start position
    # for + strand reads, adjust soft clipping at the beginning of CIGAR string
    if strand:
        result = re.match('^([0-9]+)S', cigar)
        if result:  
            adjust = 0 - int(result[1])
    # for - strand reads, sum M/D/N/S values and add that to start_pos, ignore any inserts and soft clipping at the beginning   
    else:
        result_S = re.search('([0-9]+)S$', cigar)
        result_MDN = re.findall('([0-9]+)[MDN]', cigar)
        if result_S:
            adjust += int(result_S[1])
        if result_MDN:
            adjust += sum(map(int, result_MDN))
        # subtract 1 to correct to actual position // technically not necessary
        adjust -= 1
    return start_pos + adjust

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed. validate_base_seq passed DNA and RNA tests")

    assert convert_phred("#") == 2, "convert_phred returns incorrect quality score value"
    assert convert_phred("5") == 20, "convert_phred returns incorrect quality score value"
    assert convert_phred("I") == 40, "convert_phred returns incorrect quality score value"
    print("Passed. convert_phred correctly converted phred scores")

    assert qual_score("FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@") == 37.62105263157895, "qual_score does not return correct average quality score value"
    print("Passed. qual_score correctly calculated  average quality score ")

    assert gc_content("ATGCGCGCTTAATTAA") == 0.375 , "gc_content does not return correct value"
    print("Passed. gc_content correctly calculated GC content")

    assert calc_N50([4,2,3,4,5,1,1]) == 4, "calc_N50 does not return correct N50 value"
    print("Passed. calc_N50 correctly calculated N50")

    assert rev_comp("ATTGGC") == "GCCAAT", "rev_comp does not return correct reverse complement string"
    print("Passed. rev_comp correctly reverse complemented string")

    assert meets_Qcutoff("II>III") == False, "meets_Qcutoff does not correctly return False"
    assert meets_Qcutoff("IIIIII") == True, "meets_Qcutoff does not correctly returns True"
    print("Passed. meets_Qcutoff correctly performed quality score cutoff checks")

    assert get_strand(16) == False, "get_strand does not return correct strand ("
    assert get_strand(0) == True, "get_strand does not return correct strand"
    print("Passed. get_strand correctly identified + and - strands")

    assert get_start(100, "10M", True) == 100, "get_start does not return correct start position"
    assert get_start(100, "2S8M", True) == 98, "get_start does not return correct start position"
    assert get_start(100, "8M2S", True) == 100, "get_start does not return correct start position"
    assert get_start(100, "10M", False) == 109, "get_start does not return correct start position"
    assert get_start(100, "2S8M", False) == 107, "get_start does not return correct start position"
    assert get_start(100, "8M2S", False) == 109, "get_start does not return correct start position"
    assert get_start(100, "5M2D5M", False) == 111, "get_start does not return correct start position"
    assert get_start(100, "5M2I5M", False) == 109, "get_start does not return correct start position"
    assert get_start(100, "5M10N5M", False) == 119, "get_start does not return correct start position"
    assert get_start(100, "2S5M2D2I2N5M2S", False) == 115, "get_start does not return correct start position"
    print("Passed. get_start correctly adjusted start positions")
