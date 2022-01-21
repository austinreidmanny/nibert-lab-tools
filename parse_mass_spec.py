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
from matplotlib import pyplot

from Bio import SeqRecord
from Bio import Seq
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline


# ---------------------------------------------------------------------------------------------------------------------
# Read inputs
# ---------------------------------------------------------------------------------------------------------------------
def read_inputs():

    # Initialize argparse
    parser = argparse.ArgumentParser()

    # Read in the peptides files
    parser.add_argument('--peptides', type=argparse.FileType('r'), required=True)

    # Read in the protein reference file
    parser.add_argument('--reference', type=argparse.FileType('r'), required=True)

    # Save the provided inputs
    args = parser.parse_args()
    return args


# ---------------------------------------------------------------------------------------------------------------------
# Extract peptide strings from peptide file
# ---------------------------------------------------------------------------------------------------------------------
def extract_peptide_strings(peptides_file):

    # Create an empty list to store all of the peptides
    cleaned_peptides = []

    # Setup an output filename to save the fasta of the cleaned peptides
    cleaned_peptides_filename = 'cleaned_peptides.fasta'

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
        cleaned_peptides.append(SeqRecord.SeqRecord(Seq.Seq(peptide), id=peptide_id))

    # Save cleaned peptides to a fasta file
    SeqIO.write(cleaned_peptides, cleaned_peptides_filename, 'fasta')

    return cleaned_peptides_filename


# ---------------------------------------------------------------------------------------------------------------------
# Take in the cleaned peptides fasta file and align peptides to each other
# ---------------------------------------------------------------------------------------------------------------------
def align_sequences(reference_filename, peptides_fasta):

    # Read in reference sequence and cleaned peptides and combine into one master fasta
    ref_protein = SeqIO.parse(reference_filename, 'fasta')
    peptides = SeqIO.parse(peptides_fasta, 'fasta')
    SeqIO.write(ref_protein, open('combined_proteins.fasta', 'w'), 'fasta')
    SeqIO.write(peptides, open('combined_proteins.fasta', 'a'), 'fasta')

    # Perform multiple sequence alignment using MAFFT
    mafft_cline = MafftCommandline(input='combined_proteins.fasta')
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), 'fasta')

    # Write alignment to a fasta file
    AlignIO.write(align, 'aligned_peptides.fasta', 'fasta')

    # Delete the concatenated fasta file, no use for it post-analysis
    os.remove('combined_proteins.fasta')

    # Return the alignment for visualization
    return align


# ---------------------------------------------------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------------------------------------------------
def visualize_alignment(multiple_sequence_alignment):

    """Iterate through the alignment, sequence by sequence
    If gap (-), that position is a 0, if aa residue, that position is a >1;
    Now, I want non-overlapping lines, such that every sequence gets its own line;
    If it's an aa residue, that position will get assigned the sequence_counter value,
    which increases by 1 every time we look at the next sequence"""

    msa = multiple_sequence_alignment

    sequence_counter = 0
    my_df = pandas.DataFrame(columns=['Sequence_ID', 'Position', 'Value'])

    for record in msa:

        sequence_counter += 1
        residue_position = 0

        seq_id = record.id
        sequence = record.seq

        for residue in sequence:
            residue_position += 1

            if residue == "-":
                residue_value = 0
            else:
                residue_value = sequence_counter

            my_df = my_df.append(
                {'Sequence_ID': seq_id,
                 'Position': residue_position,
                 'Value': residue_value},
                ignore_index=True
            )

    my_df.plot(kind='scatter', x='Position', y='Value')
    pyplot.show()


# ---------------------------------------------------------------------------------------------------------------------
# Main code
# ---------------------------------------------------------------------------------------------------------------------
def main():

    # Take in the user-provided inputs
    args = read_inputs()
    peptides_file = args.peptides
    reference_file = args.reference

    # Run the analysis
    cleaned_peptides = extract_peptide_strings(peptides_file)
    multiple_sequence_alignment = align_sequences(reference_file, cleaned_peptides)
    visualize_alignment(multiple_sequence_alignment)


# ---------------------------------------------------------------------------------------------------------------------
# Execute code
# ---------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
