#!/usr/bin/env python3

# -------------------------------------------------------------------------------------------------------------------- #
# sequence_length.py
# -------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------------------------------- #
# Objective:
#   This file will count the number of nucleotides in a fasta file and return the length of the sequence to the user
# -------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------------------------------- #
# Usage:
# ./sequence_length.py input.fasta
# -------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------------------------------- #
# Import necessary libraries
# -------------------------------------------------------------------------------------------------------------------- #
import sys

# -------------------------------------------------------------------------------------------------------------------- #
# Create all functions
# -------------------------------------------------------------------------------------------------------------------- #
def read_fasta_file():
    """
    This function will read in the fasta file provided by the user
    :param fasta_file: a text file containing a single nucleotide sequence with an optional header on its own line (beginning with '>')
    :return: a handle that streams the file
    """

    # Check to make sure the user has provided an input sequence in the form of a fasta file
    if len(sys.argv) < 1:
        print("Usage: sequence_length.py input.fasta")
        print("Exiting with error 1...")
        exit(1)

    fasta_file = sys.argv[1]
    fasta_handle = open(fasta_file, "r")
    return fasta_handle


def clean_input_sequence(fasta_handle):
    """
    This function will take in a data handle containing the contents of the input file and will return a cleaned
    sequence containing just nucleotide residues with no header or whitespace

    :param fasta_handle: data handle containing a nucleotide sequence, optionally with a header, line breaks are allowable
    :return: a string of nucleotide residues
    """

    # Read in fasta handle line by line, stripping whitespace and discarding any header line
    nt_sequence = ""
    for line in fasta_handle:
        line = line.strip()
        if line.startswith(">"):
            continue
        nt_sequence = nt_sequence + line

    return nt_sequence


def determine_length(nt_sequence):
    """
    This function will take in a nucleotide string, calculate its length, and return the integer length of nucleotides

    :param nt_sequence: a string of nucleotides
    :return: an integer of nucleotide length
    """

    nt_length = len(nt_sequence)
    return nt_length


def main():
    """
    With all functions written, call them each step-wise, and print to the user a simple nucleotide length

    :return: a call to print the nucleotide length
    """

    fasta_handle = read_fasta_file()
    nt_sequence = clean_input_sequence(fasta_handle)
    nt_length = determine_length(nt_sequence)
    print("Length of input sequence: {} nucleotides".format(nt_length))


# -------------------------------------------------------------------------------------------------------------------- #
# Run program
# -------------------------------------------------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()

