#!/usr/bin/env python3
# coding: utf-8
"""
Extract the best scoring BUSCO for the provided busco id
"""
"""
    !!! Since the original script was written for output format of BUSCO V3 (see https://busco.ezlab.org/busco_userguide.html#:~:text=portion%20of%20the-,BUSCO%C2%A0v3,-(PMID%3A%2029220515)%20paper ).
    !!! This script has been modfied for custom output of BUSCO v5 at 2022-04-10.
"""
import os

def fetch(busco_id, mode, folder, run_name, gene_set=None):
    """
    :param busco_id: the busco id
    :param mode: the mode, tran or prot
    :param folder: the run folder
    :param run_name: the run name
    :param gene_set: the path to the gene set fasta file
    :return: the sequence
    """

    if folder[-1] != '/':
        folder += '/'

    # load the full table
    full_table = {}
    full_table_file = open('%s/run_%s/full_table.tsv' %
                           (folder, run_name), 'r')
    for line in full_table_file:
        if not line.startswith('#'):
            if line.strip().split()[0] not in full_table:
                full_table.update(
                    {line.strip().split()[0]: [line.strip().split()[1:]]})
            else:
                full_table[line.strip().split()[0]].append(
                    line.strip().split()[1:])

    if mode == 'prot':
        return _fetch_prot(busco_id, folder, full_table, gene_set)
    elif mode == 'tran':
        return _fetch_tran(busco_id, folder, full_table)
    else:
        print('Wrong mode specified')


def _fetch_best(busco_id, folder, full_table,run_name):
    """
    :return: the sequence id to retrieve, with the best score
    """
    result_scores = []
    good_results = []
    for entry in full_table[busco_id]:
        if len(entry)>=2:
            result_scores.append(float(entry[2]))
    # open all hmm result file for this busco id
    for file in os.listdir('%s/run_%s/hmmer_output/initial_run_results' % (folder, run_name)):
        if file.startswith(busco_id):
            for line in open('%s/run_%s/hmmer_output/initial_run_results/%s' % (folder, run_name, file), 'r'):
                if not line.startswith('#') and len(line.strip()) > 1:
                    if float(line.split()[7]) in result_scores:
                        good_results.append(line)

    # Recheck that the best result has the longest protein, should be the best score as well. 
    # Edit: in fact no, for transcriptomes it is not since the 6 frames translation are evaluated and only one is correct, so the warning below are not really useful, don't focus too much on it.
    id_to_return = None
    best_length = 0
    best_score = 0
    result_with_best_length = None
    result_with_best_score = None

    for result in good_results:
        if float(result.split()[7]) > best_score:
            id_to_return = result.split()[0]
            best_score = float(result.split()[7])
            result_with_best_score = result.split()[0]
        if int(result.split()[16]) - int(result.split()[15]) > best_length:
            best_length = int(result.split()[16]) - int(result.split()[15])
            result_with_best_length = result.split()[0]

    #if result_with_best_length != result_with_best_score:
    #    _logger.warning('best score and best length are not the same for %s in %s' % (busco_id, folder))
    #else:
    #    _logger.debug('best score and best length are the same for %s in %s' % (busco_id, folder))

    return id_to_return


def _fetch_tran(busco_id, folder, full_table, run_name):
    """
    :param busco_id:
    :param folder:
    :param full_table:
    :return: the best sequence for a BUSCO id, extracted from a transcriptome BUSCO run
    """

    id_to_return = _fetch_best(busco_id, folder, full_table, run_name)

    # Now return the correct sequence
    sequences = open('%s/run_%s/translated_proteins/%s.faa' %
                     (folder, run_name, '_'.join(id_to_return.split('_')[:-1])), 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip() == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip()
    return sequence


def _fetch_prot(busco_id, folder, full_table, gene_set, run_name):
    """
    :param busco_id:
    :param folder:
    :param full_table:
    :gene_set:
    :return: the best sequence for a BUSCO id, extracted from a protein BUSCO run     """
    id_to_return = _fetch_best(busco_id, folder, full_table, run_name)
    # Now return the correct sequence
    sequences = open('%s' % gene_set, 'r')
    sequence = ''
    found = False
    for line in sequences:
        if found and line.startswith('>'):
            break
        elif line.strip().split(' ')[0] == '>%s' % id_to_return:
            found = True
            sequence += '%s\n' % line.strip()
        elif found:
            sequence += '%s' % line.strip()
    return sequence
