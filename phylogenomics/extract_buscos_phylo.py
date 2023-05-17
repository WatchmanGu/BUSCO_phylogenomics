#!/usr/bin/env python3
# coding: utf-8
"""
.. module::
   :synopsis: Extract the best BUSCOs sequence for which a single copy complete sequence is found in all species
.. moduleauthor:: Mathieu Seppey <mathieu.seppey@unige.ch>
.. versionadded:: 1.0

From BUSCO run folders, this code extracts complete single copy BUSCO genes

You will need to do something or edit the script where you find !!! >

"""
"""
    !!! Since the original script was written for output format of BUSCO V3 (see https://busco.ezlab.org/busco_userguide.html#:~:text=portion%20of%20the-,BUSCO%C2%A0v3,-(PMID%3A%2029220515)%20paper ).
    !!! This script has been modfied for custom output of BUSCO v5 at 2022-04-10.
"""

import os
import fetch_best_sequence
import re

# Config here
# !!! > you need to untar unzip the files in the run folder before running this script, if you used the -z option of BUSCO
# - single_copy_busco_sequences
# - augustus_output/extracted_proteins

wd = '/Path/to/your/BUSCO/protein_busco'
# !!! > output of this script will be here. Create the folder before running the script.
output = '/Path/to/your/BUSCO/protein_busco/extracted'
# !!! > Protein set for each species (which is the input file for BUSCO), in fasta, should be placed here
genesets = '/Path/to/your/BUSCO/pep'

if wd[-1] != '/':
    wd += '/'

if output[-1] != '/':
    output += '/'

if genesets[-1] != '/':
    genesets += '/'

def mapping_name_get(mapping_name_file):
    # !!! > Write the following mapping name file manually.
    # key = what is the name of your BUSCO result (for each sample) folders.
    # 1st entry, the (sample) name you want to see in the tree
    # If protein, write the name of the fasta containing the proteome for the sample as 3rd entry (就是每个物种用于BUSCO分析的蛋白组序列文件).
    # If transcriptome, just write tran, no 3rd entry, the sequence exists in the busco run folder
    # run_for_phylogeny_TEST => 'for_phylogeny_TEST': ['TEST', 'prot', 'protein_file_for_TEST.faa']
    # mapping_name = {
    # '': ['', 'prot',''],
    # '': ['', 'prot',''],
    # '': ['', 'prot',''],
    # '': ['', 'tran'],
    # '': ['', 'tran'],
    # '': ['', 'tran']
    # }
    # mapping_name = {
    # 'for_phylogeny_TEST': ['TEST', 'prot', 'TEST.faa'],
    # 'for_phylogeny_TEST_tran': ['TEST_tran', 'tran'],
    # }
    mapping_name_dict = {}
    with open(mapping_name_file, "r") as infile:
        for line in infile:
            line_fields = line.strip().split("\t")
            mapping_name_dict[line_fields[0]] = line_fields[1:]
    return mapping_name_dict

mapping_name = mapping_name_get(
    # !!! > Remember to change the path to your mapping_name_file.txt
    "/Path/to/your/mapping_name_file.txt")

# First, clean the output files.
for file in os.listdir(output):
    if '.faa' in file:
        os.remove('%s%s' % (output, file))

species = {}

# Then init a dict for each species and each kind of prediction (C,D,F,M)
for content in os.listdir(wd):
    if content in list(mapping_name.keys()):
        species.update(
            {content: {'Complete': [], 'Duplicated': [], 'Fragmented': [], 'Missing': []}})

# Determine which BUSCO to keep
complete = set([])
for a_species in species:
    try:
        full_table = open('%s%s/run_liliopsida_odb10/full_table.tsv' % (wd, a_species), 'r')
    except FileNotFoundError as e:
        continue
    complete_tmp = []
    for full_table_line in full_table:
        try:
            # !!! >  if you don't want the best duplicated
            if 'Complete' in full_table_line.split()[1]:
                # if 'Complete' in full_table_line.split()[1] or 'Duplicated' in full_table_line.split()[1]: # if you want to include the best scoring duplicated. In our experience, it is less reliable, but we have not really benchmarked that.
                complete_tmp.append(full_table_line.split()[0])
        except IndexError:
            pass
    if complete:
        # complete = complete.intersection(set(complete_tmp))
        complete = complete.union(set(complete_tmp))
    else:
        complete = set(complete_tmp)
# print(complete)
# Now using this list, retrieve all sequences, one file for each BUSCO. When duplicate, keep the best score
for busco in complete:
    output_file = open('%s%s.faa' % (output, busco), 'a')
    for species in mapping_name:
        if mapping_name[species][1] == 'prot':
            seq = fetch_best_sequence.fetch(busco, 'prot', '%s%s' % (wd, species),
                                            'liliopsida_odb10', # !!! > depaned on your BUSCO database
                                            '%s%s' % (genesets, mapping_name[species][2]))
        else:
            seq = fetch_best_sequence.fetch(busco, 'tran', '%s%s' % (wd, species),
                                            'liliopsida_odb10')

        search_string = r'^>.*\n'
        seq = re.sub(search_string, '>%s\n' % mapping_name[species][0], seq)
        if len(seq) >2:
            output_file.write('%s\n' % seq)
