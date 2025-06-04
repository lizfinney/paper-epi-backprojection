#!/usr/bin/env python
# coding: utf-8
# %%

import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance
import os
import subprocess
import shutil
import data_processing as dp
import pandas as pd

def find_site_index_file(filepath):
    """ Given a sequence file find the correct corresponding file that has the site names"""
    directory, file = os.path.split(filepath)
    if file.find('---')==-1:
        return filepath[:-4] + '-sites.csv'
    else:
        return filepath[:filepath.find('---')] + '-sites.csv'
    
def get_data(file, get_seqs=True):
    """ Given a sequence file, get the sequences, dates, and mutant site labels"""
    print(file)
    data      = pd.read_csv(file)
    if not get_seqs:
        sequences = None
    else:
        sequences  = np.array([np.array(list(i)) for i in list(data['sequence'])])
    dates     = list(data['date'])
    sub_dates = list(data['submission date'])
    accessions = list(data['accession'])
    index_file = find_site_index_file(file)
    index_data = pd.read_csv(index_file)
    mut_sites  = list(index_data['mutant_sites'])
    ref_sites  = list(index_data['ref_sites'])
    dic = {
        'ref_sites' : ref_sites,
        'mutant_sites' : mut_sites,
        'sequences' : sequences,
        'times' : dates,
        'submission_dates' : sub_dates,
        'accessions' : accessions
    }
    return dic


def impute_gaps(msa, counts, gaps, ref_idxs, MAX_GAP_FREQ):
    """ Impute ambiguous nucleotides with the most frequently observed ones in the alignment. """
    
    # Find site numbers that don't have known deletions
    mask     = ~np.isin(ref_idxs, gaps)
    col_idxs = np.arange(len(ref_idxs))[mask]
    
    print(np.arange(len(ref_idxs))[np.isin(ref_idxs, gaps)])
    #print('%d ambiguous nucleotides across %d sequences' % (np.sum(ambiguous), len(msa)))
    idxs_drop = []
    for i in col_idxs:
        if counts[i][0] / len(msa) > MAX_GAP_FREQ:
            idxs_drop.append(i)
    
    for i in range(len(msa)):
        gap_sites = np.array(list(msa[i]))=='0'
        idxs      = np.arange(len(msa[i]))[gap_sites]
        idxs      = [k for k in idxs if (k in col_idxs and k not in idxs_drop)]
        for j in idxs:
            new = np.argmax(counts[j, 1:]) + 1
            msa[i][j] = new
            
    idxs_keep = np.array([i for i in np.arange(len(ref_idxs)) if i not in idxs_drop])
    print(idxs_drop)
    print(idxs_keep)
    ref_idxs  = np.array(ref_idxs)[idxs_keep]
    new_msa   = []
    for seq in msa:
        new_msa.append(np.array(seq)[idxs_keep])

    return new_msa, ref_idxs


def main(args):
    """ Impute gaps """
    
    parser = argparse.ArgumentParser(description='Selection coefficients inference')
    parser.add_argument('-o',             type=str,    default=None,           help='output directory')
    parser.add_argument('--mut_dir',      type=str,    default=None,           help='directory containing mutant counts for each site')
    parser.add_argument('--input',        type=str,    default=None,           help='input file containing sequence data')
    parser.add_argument('--refFile',      type=str,    default='ref-index.csv',help='the file containing the reference index and nucleotides')
    parser.add_argument('--gapFile',      type=str,    default='gap-list.npy', help='the file containing sites that have known deletions')
    parser.add_argument('--max_gap_freq', type=float,  default=0.95,           help='the maximum frequency of gaps allowed before dropping a site')
    
    arg_list = parser.parse_args(args)
    
    out_dir   = arg_list.o        
    mut_dir   = arg_list.mut_dir
    file      = arg_list.input
    max_gap_freq = arg_list.max_gap_freq
    
    ref_df = pd.read_csv(arg_list.refFile)
    index  = list(ref_df['ref_index'])
    nucs   = np.array(ref_df['nucleotide'])
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    gap_list = np.load(arg_list.gapFile, allow_pickle=True)

    file_tail = os.path.split(file)[-1]
    print(file)

    # load mutant count file
    mut_count_data  = np.load(os.path.join(mut_dir, file_tail), allow_pickle=True)
    mut_counts = mut_count_data['mut_counts']
    seq_data   = get_data(file)
    refs       = np.array(seq_data['ref_sites'])
    muts       = np.array(seq_data['mutant_sites'])
    times      = seq_data['times']
    seq        = seq_data['sequences']
    seqs_gap   = np.array(seq_data['sequences'])


    seqs_impute, refs_impute = impute_gaps(seqs_gap, mut_counts, gap_list, refs, max_gap_freq)

    seqs_new = [''.join(list(i)) for i in seqs_impute]
  
    out_file = os.path.join(out_dir, file_tail)
    
    data = {
        'accession' : seq_data['accessions'],
        'submission_date' : seq_data['submission_dates'],
        'date' : times,
        'sequence' : seqs_new
    }
    df = pd.DataFrame(data=data)
    df.to_csv(out_file, index=False)
        
    sites_data = {'ref_sites' : refs_impute}
    df2 = pd.DataFrame(data=sites_data)
    df2.to_csv(out_file[:-4] + '-sites.csv', index=False)
    
    
if __name__ == '__main__': 
    main(sys.argv[1:])