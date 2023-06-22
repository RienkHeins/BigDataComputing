#!/bin/bash

# create files variable
export FILES=/data/dataprocessing/rnaseq_data/Brazil_Brain/

# Using hosts.txt with assemblix 2012 only as the Brazil_Brain folder couldn't be found in other assemblix directories
# or nuc systems, bin systems wouldn't connect through ssh.

# run parralel command
find $FILES -name '*.fastq' | parallel -j 2 --sshloginfile hosts.txt --results ./outdir/ "python3 assignment4.py {}"
