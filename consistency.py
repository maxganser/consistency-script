#!/usr/bin/env python3

# The MIT License (MIT)
# Copyright (c) 2021 Thomas Huetter.
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
Program: Identify consensus and non-consensus molecular characters between 
multiple alignment files.
'''

import argparse
import csv
import os
import numpy as np
from Bio import AlignIO
import DeSignate.webapp.designate as DS

def create_position_alignment(alignment):
    """
        Given an alignment with multiple sequences, assign ascending numbers 
        representing their positions. Positions with gaps are skipped and 
        assigned 0.
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

def get_signature_characters(alignment, query_group, reference_group, k):
    """
        Calls DeSignate to detect signature characters and filters for binary 
        and asymmetric characters. Returns a list of the according positions 
        in the alignment.
    """

    signature_characters = [r[0] for r in DS.signature_character_detection(
        alignment, None, query_group, reference_group, k, None, False, True, 
        False) if r[4] == 'binary' or r[4] == 'asymmetric']

    return signature_characters


def main():
    # Parse input arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignments", type=str, nargs="+", help="<Required> \
        Path to input files that contain the alignments. The first one \
        is considered to be the reference alignment.", required=True)
    parser.add_argument("--alignment_format", type=str, default="fasta",
        help="Specify the alignment file format. Default is set to 'fasta'.")
    parser.add_argument("--query_group", type=str, help="<Required> \
        Path to a csv file that contains the species of the query group.", 
        required=True)
    parser.add_argument("--reference_group", type=str, help="<Required> \
        Path to a csv file that contains the species of the reference group.", 
        required=True)
    parser.add_argument("--k_window", type=int, default=1,
        help="Combine two noisy positions in order to find more asymmetric \
        pairs resulting in a higher number of unique characteristics. The \
        required parameter k is used for the so-called k-window. Only noisy \
        positions in a range of k are considered to be combined for analysis.")
    args = parser.parse_args()

    # List of all input alignments.
    alignments = []
    # List of all input alignments with position numbers instead of nucleotides.
    position_alignments = []
    # List of all signature characters detected by DeSignate for each alignment.
    signature_characters = []
    # List of the species in query and reference groups. Read given CSV files 
    # as input.
    query_group = []
    with open(args.query_group, newline='') as f:
        reader = csv.reader(f)
        query_group = list(reader)[0]

    reference_group = []
    with open(args.reference_group, newline='') as f:
        reader = csv.reader(f)
        reference_group = list(reader)[0]

    # Parse the alignments, translate to position alignments, and analyze 
    # their signature characters using DeSignate.
    for alignment in args.alignments:
        alignments.append(AlignIO.read(alignment, args.alignment_format))
        position_alignments.append(
            create_position_alignment(alignments[len(alignments)-1]))
        # Perform DeSignate analysis.
        signature_characters.append(get_signature_characters(
            alignments[len(alignments)-1], query_group, reference_group, 
            args.k_window))


    ############################################################################
    ####### TODO BELOW #########################################################
    ############################################################################

    # Turn alignment_vector into array
    # Columns can be selected in the arrays
    align_array_A1 = np.array(position_alignments[0], np.dtype(int))
    # print("Array shape %i by %i" % align_array_A1.shape)

    align_array_A2 = np.array(position_alignments[1], np.dtype(int))
    #print("Array shape %i by %i" % align_array_muscle.shape)

    align_array_A3 = np.array(position_alignments[2], np.dtype(int))
    #print("Array shape %i by %i" % align_array_tcoffee.shape)

    # Selected columns get new indices starting from 0 so pos-1 is needed
    # A1
    print('\nReference alignment (A1)')
    # Select positions in array and write to new list of tuples (hashable)
    sig_pos_A1 = []
    for pos in signature_characters[0]:
        sig_pos_A1.append(tuple(align_array_A1[:, (int(pos)-1)]))
    print('Number of signatures selected:', len(sig_pos_A1))

    # A2
    print('\nAlternative alignment (A2)')
    # Select positions in array and write to new list of tuples (hashable)   
    sig_pos_A2 = []
    for pos in signature_characters[1]:
        sig_pos_A2.append(tuple(align_array_A2[:, (int(pos)-1)]))
    print('Number of signatures selected:', len(sig_pos_A2))

    # A3
    print('\nAlternative alignment (A3)')
    # Select positions in array and write to new list of tuples (hashable)   
    sig_pos_A3 = []
    for pos in signature_characters[2]:
        sig_pos_A3.append(tuple(align_array_A3[:, (int(pos)-1)]))
    print('Number of signatures selected:', len(sig_pos_A3))

    # Create dictionaries to link signature positions with nucleotide positions and states
    sig_nucpos_A1 = dict(zip(signature_characters[0], sig_pos_A1))
    sig_nucpos_A2 = dict(zip(signature_characters[1], sig_pos_A2))
    sig_nucpos_A3 = dict(zip(signature_characters[2], sig_pos_A3))

    # Convert dictionaries to also get the corresponding signature positions 
    # from each alignment
    A2_dict = dict((y,x) for x,y in sig_nucpos_A2.items())
    A3_dict = dict((y,x) for x,y in sig_nucpos_A3.items())

    # Get consensus and non-consensus signatures, save to list, and print
    # Which consensus signature positions correspond to each other in each alignment? 
    consensus_signatures = []
    for k in sig_nucpos_A1:
        if sig_nucpos_A1[k] in A2_dict and sig_nucpos_A1[k] in A3_dict:
            consensus_signatures.append((int(k), 
                                              int(A2_dict[sig_nucpos_A1[k]]), 
                                              int(A3_dict[sig_nucpos_A1[k]])))
            
    non_consensus_signatures = []
    for k in sig_nucpos_A1:
        if sig_nucpos_A1[k] not in A2_dict or sig_nucpos_A1[k] not in A3_dict:
            non_consensus_signatures.append((int(k)))

    # Print results
    print('\nCONSENSUS SIGNATURE CHARACTERS')
    print('Quantity:',
          len(consensus_signatures), 
          '\nPositions:')
    print("Alignment 1  Alignment 2  Alignment 3")
    for pos1,pos2,pos3 in consensus_signatures:
        print("{:<13}{:<13}{:<13}".format(pos1,pos2,pos3))
    print ('\nNON-CONSENSUS SIGNATURE CHARACTERS',
           '\nQuantity:',
           len(non_consensus_signatures), 
          '\nPositions in reference alignment (A1):')
    for pos1 in non_consensus_signatures:
        print(pos1)

    # print(query_group)
    # print(reference_group)
    # print(alignments)
    # print(position_alignments[0][0])
    # print(signature_characters[0])


if __name__ == "__main__":
    main()
