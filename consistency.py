#!/usr/bin/env python3

# The MIT License (MIT)
# Copyright (c) 2021 Thomas Huetter & Maximilian Ganser.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

'''
Program: Identify consensus and non-consensus signature characters between 
a reference alignment and alternative alignments.
'''

import argparse
import csv
import os
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
import DeSignate.webapp.designate as DS
import pandas as pd

def create_position_alignment(alignment):
    """
        Given an alignment with multiple sequences, assign ascending numbers 
        to the nucleotides of each sequence representing their positions. 
        Positions with gaps are skipped and assigned 0.
    """

    position_alignment = []
    for record in alignment:
        position_sequence = []
        pos = 1
        for c in record:
            if c == '-':
                # Assign 0 and skip gaps.
                position_sequence.append(0)
            else:
                # Assign ascending numbers.
                position_sequence.append(pos)
                pos += 1
        position_alignment.append(position_sequence)

    return position_alignment

def get_sig_chars(alignment, query_group, reference_group, k, gaps):
    """
        Calls DeSignate to detect signature characters and filters for binary 
        and asymmetric signatures. Returns a list of the corresponding 
        alignment positions.
    """

    sig_chars = [r[0] for r in DS.signature_character_detection(alignment, 
        None, query_group, reference_group, k, "signature-characters.csv", 
        True, True, gaps) if r[4] == 'binary' or r[4] == 'asymmetric']

    return sig_chars


def column(matrix, col):
    """
        Returns a given column of a two-dimensional list.
    """

    return [row[col] for row in matrix]

def print_matrix(matrix):
    """
        Prints a given 2D list (matrix) nicely.
    """

    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))

    return

def main():
    # Parse input arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignments", type=str, nargs="+", help="<Required> \
        Path to input files that contain the alignments. The first one \
        is considered to be the reference alignment.", required=True)
    parser.add_argument("--alignment_format", type=str, default="fasta",
        help="Specify the alignment file format. Default is set to 'fasta'.")
    parser.add_argument("--query_group", type=str, help="<Required> \
        Path to a csv file that contains the sequence labels of the query \
        group.", required=True)
    parser.add_argument("--reference_group", type=str, help="<Required> \
        Path to a csv file that contains the sequence labels of the reference \
        group.", required=True)
    parser.add_argument("--k_window", type=int, default=1,
        help="Combine two noisy positions in order to find more asymmetric \
        pairs resulting in a higher number of unique characteristics. The \
        required parameter k is used for the so-called k-window. Only noisy \
        positions in a range of k are considered to be combined for analysis.")
    parser.add_argument("--consider_gaps", default=True, action="store_true", 
        help="Classify signature characters with gap state.")
    args = parser.parse_args()

    # List of all input alignments.
    alignments = []
    # List of all input alignments with position numbers instead of nucleotides.
    numbered_alignments = []
    # List of all signature characters detected by DeSignate for each alignment.
    sig_chars = []
    # Stores the number of each column that belongs to a signature character.
    numbered_sig_cols = []
    # Dictionary from signature character positions to numbered column arrays.
    sig_to_cols = []
    # Dictionary from numbered column arrays to signature character positions.
    cols_to_sig = []
    # Stores consensus and non-consensus signatures.
    consensus_signatures = []
    non_consensus_signatures = []

    # List of the species in query and reference groups. Read given CSV files 
    # as input.
    query_group = []
    with open(args.query_group, newline='') as f:
        reader = csv.reader(f, skipinitialspace=True)
        query_group = list(reader)[0]

    reference_group = []
    with open(args.reference_group, newline='') as f:
        reader = csv.reader(f, skipinitialspace=True)
        reference_group = list(reader)[0]

    # Parse the alignments, translate to nucleotide position alignments, and  
    # analyze their signature characters using DeSignate.
    for alignment in args.alignments:
        # Parse input alignments and keep all sequences that are either part 
        # of the query or the reference group.
        sequences = []
        for seq in SeqIO.parse(alignment, args.alignment_format):
            if seq.name in query_group or seq.name in reference_group:
                sequences.append(seq)
        alignments.append(MultipleSeqAlignment(sequences))

        # Build an alignment that stores ascending numbers rather than 
        # characters.
        numbered_alignments.append(
            create_position_alignment(alignments[-1]))

        # Perform DeSignate analysis.
        sig_chars.append(get_sig_chars(
            alignments[-1], query_group, reference_group, 
            args.k_window, args.consider_gaps))

        # Extract all numbers at signature character columns. Tuples are needed 
        # to hash them later on.
        numbered_sig_cols.append([])
        for pos in sig_chars[-1]:
            numbered_sig_cols[-1].append(
                    tuple(sorted(column(numbered_alignments[-1], pos-1))))

        # Link signature character positions to numbered columns.
        sig_to_cols.append(dict(zip(sig_chars[-1], numbered_sig_cols[-1])))
        # Link numbered columns to signature character positions.
        cols_to_sig.append(dict(zip(numbered_sig_cols[-1], sig_chars[-1])))

    # Get signature characters of reference alignment and save as csv file
    designate_results = [r[0] for r in DS.signature_character_detection(
        alignments[0], None, query_group, reference_group, 1,
        "signature-characters.csv", True, False, True)]

    # Calculate Shannon-entropy values of reference alignment
    shannon_entropy = [r[0] for r in DS.shannon_entropy_analysis(alignments[0],
        None, None, "shannon-entropy.csv")]
        
    # Identify consensus and non-consensus signature characters. A consensus 
    # signature character has identical column numbers to the reference 
    # alignment (i.e., index 0).
    for k, v in sig_to_cols[0].items():
        # Check whether a reference alignment signature occurs in all other 
        # alignments.
        consensus = True
        # Start at second alignment since the reference alignment is not 
        # considered.
        for cst in cols_to_sig[1:]:
            if v not in cst:
                consensus = False
        if consensus is True:
            cons_sig = [k]
            # Start at second alignment since the reference alignment is not 
            # considered.
            for cst in cols_to_sig[1:]:
                cons_sig.append(cst[v])
            consensus_signatures.append(cons_sig)
        else:
            non_consensus_signatures.append([k])


    # Save results.
    # Build header for consensus matrix.
    consensus_matrix_header = []
    for i in range(len(cols_to_sig)):
        consensus_matrix_header.append("Alignment " + str(i+1))

    writer = csv.writer(open("consensus-positions.csv", "w"))
    writer.writerow(consensus_matrix_header)
    writer.writerows(consensus_signatures)

    writer = csv.writer(open("non-consensus-positions.csv", "w"))
    writer.writerow(["Reference Alignment"])
    writer.writerows(non_consensus_signatures)

    # Build final results file
    a = pd.read_csv("signature-characters.csv")
    b = pd.read_csv("shannon-entropy.csv")
    c = pd.read_csv("consensus-positions.csv")
    merged_a = a.merge(b, on='position')
    merged_a.to_csv("designate-results.csv", index=False)
    d = pd.read_csv("designate-results.csv")
    merged_c = c.merge(d, left_on='Alignment 1', right_on='position')
    merged_c.to_csv("consensus-sigchars.csv", index=False)
    
    # Print results
    print('CONSENSUS SIGNATURE CHARACTERS')
    print('Quantity:' + str(len(consensus_signatures)))
    print('Positions:')
    print_matrix(consensus_signatures)
    print()
    print('NON-CONSENSUS SIGNATURE CHARACTERS')
    print('Quantity:' + str(len(non_consensus_signatures)))
    print('Positions in the reference alignment:')
    for pos in non_consensus_signatures:
        print(pos[0])

if __name__ == "__main__":
    main()
