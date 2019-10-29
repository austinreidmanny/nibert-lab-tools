#!/bin/bash

# Launch conda environment to load tools
eval "$(conda shell.bash hook)"
conda activate env_bbsplit

# Initialize list of taxonomy IDs
if [[ "$#" -ne 1 ]]; then
    echo -e "ERROR: No list of taxonomy IDs detected. \n" \
            "Usage: $0 taxonomy_list.txt \n" \
            "EXITING..." && exit 1
fi

taxonomy_IDs=$1

# Run the BBtools taxonomy tool
taxonomy.sh \
    tree="~/Documents/Research/Tools/taxonomy/tree.taxtree.gz" \
    $(cat ${taxonomy_IDs})

