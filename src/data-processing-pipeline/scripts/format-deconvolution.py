#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance
import os
import subprocess
import shutil
#import data_processing as dp
import pandas as pd

REF_TAG = 'EPI_ISL_402125'

def smoothing_to_MPL(smooth_dir, nonunique_dir, prob_dist, nosmooth_file, out_file, timed=1):

    """ Takes the output directory containing files that have times and distributions resulting from the deconvolution method,
    and outputs a file that can be read by the inference script.
    prob_dist is a list containing the probability distribution for a single sequence.
    nosmooth file is the file containing the data before smoothing."""
    
    if timed>0:
        t_start = timer()
     
    # get the mutant sites from the non-smoothed data file

    muts = list(pd.read_csv(nosmooth_file)['ref_sites'])
    
    # load the file that identifies sequences with an index

    df_idx      = pd.read_csv(os.path.join(smooth_dir, 'sequence-index.csv'))
    idxs        = list(df_idx['index'])
    seqs_mult   = list(df_idx['sequences'])
    
    # load the file containing sequences that only appear once

    df_ind      = pd.read_csv(os.path.join(smooth_dir, 'unique-sequences.csv'))

    seqs_ind = list(df_ind['sequences'])
    times_ind   = list(df_ind['times'])
    times_ind   = [np.arange(times_ind[i] - len(prob_dist) + 1, times_ind[i] + 1) for i in range(len(times_ind))]
    
    # load the data for the sequences that appear multiple times

    times_mult  = np.empty(len(idxs), dtype=object)
    counts_mult = np.empty(len(idxs), dtype=object)
    seqs_new    = np.empty(len(idxs), dtype=object)

    # fix this line

    for file in os.listdir(nonunique_dir):

        if file.endswith('_rates.csv'):
            identifier_str = file.split('_')[0]  # Extract the part before '_rates.csv'
            identifier = int(identifier_str)
            df_temp     = pd.read_csv(os.path.join(nonunique_dir, file))
            times_temp  = list(df_temp['times'])
            counts_temp = list(df_temp['rates'])
            seq_idx     = idxs.index(identifier)
            times_mult[seq_idx]  = times_temp
            counts_mult[seq_idx] = counts_temp
    
    print(times_ind)
    print(times_mult)

    # find all times and set up arrays to be filled with sequences and counts

    times_all = np.sort(np.unique([times_ind[i][j]  for i in range(len(times_ind))  for j in range(len(times_ind[i]))] + 
                                  [times_mult[i][j] for i in range(len(times_mult)) for j in range(len(times_mult[i]))]))
    times_all = list(times_all)
    
    if timed>0:
        t_end = timer()
        print('total time = ', t_end - t_start)

    # make .csv file

    headers = 'sequence,' + ','.join(np.array(times_all, dtype=str)) + '\n'
    g = open(out_file, mode='w')
    g.write(headers)
    for i in range(len(seqs_ind)):
        time_idxs   = [times_all.index(j) for j in times_ind[i]]
        counts_temp = np.zeros(len(times_all))
        for j in range(len(time_idxs)):
            counts_temp[time_idxs[j]] += prob_dist[j]
        counts_temp = np.array(counts_temp, dtype=str)
        counts_str = ','.join(counts_temp)
        g.write(f'{seqs_ind[i]},{counts_str}\n')
    for i in range(len(seqs_mult)):
        time_idxs   = [times_all.index(j) for j in times_mult[i]]
        counts_temp = np.zeros(len(times_all))
        for j in range(len(time_idxs)):
            counts_temp[time_idxs[j]] += counts_mult[i][j]
        counts_temp = np.array(counts_temp, dtype=str)
        counts_str = ','.join(counts_temp)
        g.write(f'{seqs_mult[i]},{counts_str}\n')
    g.close()
    

def main(args):
    """ Formats the regional sequence data so that it can be analyzed by the smoothing deconvolution method. """
    
    parser = argparse.ArgumentParser(description='Selection coefficients inference')
    parser.add_argument('-o',             type=str,    default=None,                               help='output directory')
    parser.add_argument('--input_dir',    type=str,    default=None,                               help='submit directory')
    parser.add_argument('--sites_dir',    type=str,    default='impute-gaps',                      help='sites file directory')
    parser.add_argument('--location',     type=str,    default=None,                               help='name of region')
    parser.add_argument('--prob_dist',    type=str,    default='single-sequence-distribution.csv', help='probability distribution for single sequence')
    
    arg_list = parser.parse_args(args)
    
    out_dir   = arg_list.o        
    input_dir = arg_list.input_dir
    sites_dir = arg_list.sites_dir
    location  = arg_list.location
    prob_dist = arg_list.prob_dist
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    smooth_dir = os.path.join(input_dir, 'genome-unique', location)
    nonunique_dir = os.path.join(smooth_dir, f'bp-output-{location}')

    prob_dist_values = pd.read_csv(prob_dist)['bp_pmf_value']

    sites_file = location + '-sites.csv'
    sites_file_path = os.path.join(sites_dir, sites_file)

    out_file = location + '-formatted.csv'
    out_file_path = os.path.join(out_dir, out_file)

    print(smooth_dir)
    print(nonunique_dir)
    print(sites_file_path)
    print(out_file_path)
    

    smoothing_to_MPL(smooth_dir, nonunique_dir, prob_dist_values, sites_file_path, out_file_path)

    
if __name__ == '__main__': 
    main(sys.argv[1:])
    
    