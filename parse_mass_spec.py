#!/usr/bin/env/python3

# =====================================================================================================================
# parse mass spec
# =====================================================================================================================

# ---------------------------------------------------------------------------------------------------------------------
# Objective
#
# This script will read in a 'peptides.csv' mass spec results file from the Harvard Medical School Mass Spec Core Lab.
# It will extract the Peptides column, process each peptide, and then align them to a reference protein.
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# Usage
#
# ./parse_mass_spec.py --peptides peptides.csv --reference protein.fasta
# ---------------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------
# Import libraries
# ---------------------------------------------------------------------------------------------------------------------
import argparse
import io
import os
import pandas

import Bio.SeqRecord
from Bio import Seq
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline

# ---------------------------------------------------------------------------------------------------------------------
# Read inputs
# ---------------------------------------------------------------------------------------------------------------------
# Initialize argparse
parser = argparse.ArgumentParser()

# Read in the peptides files
parser.add_argument('--peptides', type=argparse.FileType('r'), required=True)

# Read in the protein reference file
parser.add_argument('--reference', type=argparse.FileType('r'), required=True)

# Save the provided inputs
args = parser.parse_args()
peptides_file = args.peptides
reference_sequence = args.reference


# ---------------------------------------------------------------------------------------------------------------------
# Extract peptide strings from peptide file
# ---------------------------------------------------------------------------------------------------------------------
def extract_peptide_strings():

    # Create an empty list to store all of the peptides
    cleaned_peptides = []

    # Read the file as a DataFrame and extract the uncleaned peptide strings
    peptides_df = pandas.read_csv(peptides_file, sep=',', usecols=['ScanF', 'Peptide']).astype('string')
    peptides_with_ids = dict(peptides_df.values)

    # Read in the peptide strings one-by-one, clean them up, and save to new dictionary
    for peptide_id in peptides_with_ids.keys():

        # Read in the entry
        peptide = peptides_with_ids[peptide_id]

        # Remove the first and last predicted residues
        peptide = peptide.split('.')[1]

        # Remove any "*" characters
        peptide = peptide.replace('*', '')

        # Save the cleaned peptide string
        cleaned_peptides.append(Bio.SeqRecord.SeqRecord(Seq.Seq(peptide), id=peptide_id))

    # Save cleaned peptides to a fasta file
    SeqIO.write(cleaned_peptides, 'cleaned_peptides.fasta', 'fasta')


# ---------------------------------------------------------------------------------------------------------------------
# Take in the cleaned peptides fasta file and align peptides to each other
# ---------------------------------------------------------------------------------------------------------------------
def align_sequences(peptides_fasta):

    # Read in reference sequence and cleaned peptides and combine into one master fasta
    ref_protein = SeqIO.parse(reference_sequence, 'fasta')
    peptides = SeqIO.parse(peptides_fasta, 'fasta')
    SeqIO.write(ref_protein, open('combined_proteins.fasta', 'w'), 'fasta')
    SeqIO.write(peptides, open('combined_proteins.fasta', 'a'), 'fasta')

    # Perform multiple sequence alignment using MAFFT
    mafft_cline = MafftCommandline(input='combined_proteins.fasta')
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), 'fasta')
    print(align)

    AlignIO.write(align, 'aligned_peptides.fasta', 'fasta')
    os.remove('combined_proteins.fasta')

# ---------------------------------------------------------------------------------------------------------------------
# Main code
# ---------------------------------------------------------------------------------------------------------------------
# Run the main function
if __name__ == '__main__':
    extract_peptide_strings()
    align_sequences('cleaned_peptides.fasta')
