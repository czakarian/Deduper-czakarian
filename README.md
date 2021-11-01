# Deduper - A Reference Based PCR Duplicate Removal Tool

Given a sorted SAM file of uniquely mapped reads, will remove all PCR duplicates and retain only a single copy of each read. 

### Input
1. SAM file sorted by chromosome/position (samtools sort)
2. Text file with a list of UMIs (will assume random UMIs if file not provided)

### Output
1. SAM file with PCR duplicates removed (filename}_deduped.sam)
2. Summary statistics (# of dups removed, unique reads per chromosome)

- Includes options for:
    - Single-end vs paired-end // not implemented yet
    - Known UMIs vs randomers (error correction?) // not implemented yet
    - Choice of duplicate written to file (first encountered or highest quality) // not implemented yet

- argparse options:
    - ```-f```, ```--file```: required arg, absolute file path of SAM file
    - ```-p```, ```--paired```: optional arg, designates file is paired end
    - ```-u```, ```--umi```: optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)
