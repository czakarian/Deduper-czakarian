# Deduper - A Reference Based PCR Duplicate Removal Tool

Given a sorted SAM file of uniquely mapped reads, this program will remove all PCR duplicates, retaining only a single copy of each read. 
PCR duplicates will have the same UMI and start position (after correcting for soft clipping) on the same chromosome and strand.

### Input
1. SAM file sorted by chromosome/position (samtools sort)
2. Text file with a list of UMIs (will assume random UMIs if file not provided)

- argparse options:
    - ```-f```, ```--file```: required arg, absolute file path of SAM file
    - ```-p```, ```--paired```: optional arg, designates file is paired end (not yet implemented)
    - ```-u```, ```--umi```: optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)
    - ```-q```, ```--quality```: optional arg, will output duplicate read of highest average quality (if unset, will output first duplicate encountered)

### Output
1. SAM file with PCR duplicates removed 
    - {filename}_deduped.sam
2. Summary statistics 
    - number of dups removed
    - number of reads with invalid UMIs 
    - number of unique reads per chromosome
