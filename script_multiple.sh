#!/bin/bash

##Note: the GTF, FASTA and FASTQ files should have the name name 
##(but not the same extension though..)
for FILE in ./ref-genomes/*.fa
do
    ## The '&' is important 'cause it allows parallelism (which is nice)
    ## Be careful that actually this works if all files are single-end 
    ## You should just modify script.py in the case of having SE and PE files
    python script $FILE &
done

## BIBLIOGRAPHY:
## To know the difference between PE (paired end) and SE (single end) reads:
# http://seqanswers.com/forums/showthread.php?t=503
# In fact, PE reads are just reads from the sequencing of both ends of the DNA molecule
# (file1 is the sequencing of the molecule when starting from the 3' end, and file2 from the 5' end)
# SE reads are obtained by sequencing from one end of the molecule
# (as far as I understand this)