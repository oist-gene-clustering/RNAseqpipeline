#!/bin/bash

##Note: the GTF, FASTA and FASTQ files should have the name name 
##(but not the same extension though..)
for FILE in ./ref-genomes/*.fa
do
    ## The '&' is important 'cause it allows parallelism (which is nice)
    ## Be careful that actually this works if all files are single-end 
    ## You should just modify script-sc.py in the case of having SE and PE files
    python script-sc $FILE &
done
