#!/bin/bash

# Usage: ./deduplicate.sh INPUT.FASTA

INPUT_FASTA=$1
OUTPUT_FILE=${INPUT_FASTA%.*}'.deduplicated.fasta'

# Check that an input file was given
if [[ "$#" -ne 1 ]] ; then
    echo -e "\n ERROR: Input fasta file was not provided. \n\n" \
            "Usage: $0 INPUT.fasta \n\n" \
            "Exiting... \n"
    exit 1
fi

# Make sure that `seqtk` is installed
command -v seqtk > /dev/null || \
{ echo -e "ERROR: `seqtk` is required. \n"
        "Please install from github.com/lh3/seqtk and try again" && \
  exit 2
}                  
        
# Take the input fasta, transform it into tab delimited format
# Sort by the first col (names), keeping only the seqs with unique names
# Then turn it back into fasta

seqtk subseq -t ${INPUT_FASTA} <(grep "^>" ${INPUT_FASTA} | \
cut -d ">" -f 2) | cut -f 1,3 | \
sort -k1,1 -u | \
awk '{print ">"$1"\n"$2}' > \
${OUTPUT_FILE}

