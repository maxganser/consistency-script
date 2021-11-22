# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 14:39:42 2021

@author: Maximilian Ganser
"""

#import argparse
import numpy as np
from Bio import AlignIO

# Parse input alignments and sort
alignment_01 = AlignIO.read("groupXY_reference_alignment.fasta", "fasta")
alignment_01.sort()
alignment_02 = AlignIO.read("groupXY_alignment_01.fasta", "fasta")
alignment_02.sort()
alignment_03 = AlignIO.read("groupXY_alignment_02.fasta", "fasta")
alignment_03.sort()

# Function to convert every sequence record in the alignment to numbers
# (gap = 0, nucleotides are numbered consecutively)
# Reference alignment (A1)
alignment_vector_01 = []
for record in alignment_01:
    seq_list_01 = []
    i = 1
    for c in record:
        if c == '-':
            seq_list_01.append('0')
        else:
            seq_list_01.append(str(i))
            i += 1
    alignment_vector_01.append(seq_list_01)

# Alternative alignment (A2)    
alignment_vector_02 = []
for record in alignment_02:
    seq_list_02 = []
    i = 1
    for c in record:
        if c == '-':
            seq_list_02.append('0')
        else:
            seq_list_02.append(str(i))
            i += 1
    alignment_vector_02.append(seq_list_02)

# Alternative alignment (A3)    
alignment_vector_03 = []
for record in alignment_03:
    seq_list_03 = []
    i = 1
    for c in record:
        if c == '-':
            seq_list_03.append('0')
        else:
            seq_list_03.append(str(i))
            i += 1
    alignment_vector_03.append(seq_list_03)

# Select signatures
# A1
f = open('reference_signature_positions.txt', 'r')
sigp_A1 = f.read().splitlines()

# A2
f = open('alignment1_signature_positions.txt', 'r')
sigp_A2 = f.read().splitlines()

# A3
f = open('alignment2_signature_positions.txt', 'r')
sigp_A3 = f.read().splitlines()

# Turn alignment_vector into array
# Columns can be selected in the arrays
align_array_A1 = np.array(alignment_vector_01, np.dtype(int))
#print("Array shape %i by %i" % align_array_A1.shape)

align_array_A2 = np.array(alignment_vector_02, np.dtype(int))
#print("Array shape %i by %i" % align_array_muscle.shape)

align_array_A3 = np.array(alignment_vector_03, np.dtype(int))
#print("Array shape %i by %i" % align_array_tcoffee.shape)

# Selected columns get new indices starting from 0 so pos-1 is needed
# A1
print('\nReference alignment (A1)')
# Select positions in array and write to new list of tuples (hashable)
sig_pos_A1 = []
for pos in sigp_A1:
    sig_pos_A1.append(tuple(align_array_A1[:, (int(pos)-1)]))
print('Number of signatures selected:', len(sig_pos_A1))

# A2
print('\nAlternative alignment (A2)')
# Select positions in array and write to new list of tuples (hashable)   
sig_pos_A2 = []
for pos in sigp_A2:
    sig_pos_A2.append(tuple(align_array_A2[:, (int(pos)-1)]))
print('Number of signatures selected:', len(sig_pos_A2))

# A3
print('\nAlternative alignment (A3)')
# Select positions in array and write to new list of tuples (hashable)   
sig_pos_A3 = []
for pos in sigp_A3:
    sig_pos_A3.append(tuple(align_array_A3[:, (int(pos)-1)]))
print('Number of signatures selected:', len(sig_pos_A3))

# Create dictionaries to link signature positions with nucleotide positions and states
sig_nucpos_A1 = dict(zip(sigp_A1, sig_pos_A1))
sig_nucpos_A2 = dict(zip(sigp_A2, sig_pos_A2))
sig_nucpos_A3 = dict(zip(sigp_A3, sig_pos_A3))

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