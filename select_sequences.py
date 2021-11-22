# -*- coding: utf-8 -*-
"""
Created on Fri Jul 30 11:19:55 2021

@author: Maximilian Ganser
"""
#
from Bio import SeqIO 
#input_file = "reference_alignment.fasta"
#input_file = "alternative_alignment_01.fasta"
input_file = "alternative_alignment_02.fasta"

id_file = "group_selection.txt"

#output_file = "groupXY_reference_alignment.fasta"
#output_file = "groupXY_alignment_01.fasta"
output_file = "groupXY_alignment_02.fasta"

with open(id_file) as id_handle:
    group_selection = set(line.rstrip("\n").split(None, 1)[0] for line in id_handle)
    print("found %i unique identifiers in %s" % (len(group_selection), id_file))
    
records = (r for r in SeqIO.parse(input_file, "fasta") if r.id in group_selection)
count = SeqIO.write(records, output_file, "fasta")
print("Saved %i records from %s to %s" % (count, input_file, output_file))
if count < len(group_selection):
    print("Warning %i IDs not found in %s" % (len(group_selection) - count, input_file))
    
print(records)





