import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance
import os
import subprocess
import shutil
import pandas as pd

REF_TAG = 'EPI_ISL_402125'
NUC     = ['-', 'A', 'C', 'G', 'T']

# _ _ _ Functions for reading data _ _ _ #

def find_site_index_file(filepath, trim_dir):
    """Given a sequence file, find the correct corresponding file that has the site names"""
    directory, file = os.path.split(filepath)
    parent_directory = os.path.dirname(directory)
    pipe_directory = os.path.join(parent_directory, trim_dir)
    if file.find('---') == -1:
        return os.path.join(pipe_directory, file[:-4] + '-sites.csv')
    else:
        return os.path.join(pipe_directory, file[:file.find('---')] + '-sites.csv')

def get_data(file, trim_dir, get_seqs=True):
    """ Given a sequence file, get the sequences, dates, and mutant site labels"""
    data = pd.read_csv(file)
    if not get_seqs:
        sequences = None
    else:
        sequences = np.array([list(i) for i in list(data['sequence'])])
    times = list(data['times'])
    counts = list(data['counts'])
    
    #parent_directory = os.path.dirname(os.path.dirname(file))  # two directories up (change later?)
    index_file = find_site_index_file(file, trim_dir)
    index_data = pd.read_csv(index_file)
    mut_sites = list(index_data['mutant_sites'])
    ref_sites = list(index_data['ref_sites'])
    dic = {
        'ref_sites': ref_sites,
        'mutant_sites': mut_sites,
        'sequences': sequences,
        'times': times,
        'counts': counts
    }
    return dic

def find_max_freq(file, trim_dir, out_dir, ref_seq=None, ref_labels=None, window=None, d=5):
    """ Finds the maximum frrequency that mutant alleles at each site reach."""
    df                = pd.read_csv(file)
    data              = get_data(file, trim_dir)
    sequences         = data['sequences']
    ref_sites         = data['ref_sites']
    mutant_sites_samp = ref_sites
    dates             = data['times']
    filename          = os.path.basename(file)
    sequences         = np.array([list(i) for i in sequences], dtype=int)
    site_idxs = [ref_labels.index(str(i)) for i in ref_sites]
    ref_poly  = np.array(ref_seq)[np.array(site_idxs)]

    unique_dates = np.unique(dates)
    num_seqs     = [np.sum(df[df['times']==i]['counts']) for i in unique_dates]
    for i in unique_dates:
        count_redo = df[df['times']==i]['counts']

    counts_tot = np.zeros((len(unique_dates), len(mutant_sites_samp), d))  
    
    for t in range(len(unique_dates)):

        idxs = np.where(dates == unique_dates[t])[0]

        seqs_t = sequences[idxs]

        count_idx = np.where(df['times'] == unique_dates[t])[0]

        count_redo = df[df['times'] == unique_dates[t]]['counts']  # Reassign count_redo based on the current unique date

        for i in range(len(mutant_sites_samp)):
 
            for j in range(len(seqs_t)):
                counts_tot[t, i, seqs_t[j][i]] += count_redo.iloc[j]  # Use count_redo.iloc[j] here

    counts = np.sum(counts_tot, axis=0)
    freq = np.array([counts_tot[i] / num_seqs[i] for i in range(len(counts_tot))])
    count_arr = np.array([counts_tot[i] for i in range(len(counts_tot))])

    mut_freqs  = np.zeros((len(unique_dates), len(mutant_sites_samp), d-1))
    mut_counts = []
    for i in range(len(mutant_sites_samp)):
        ref_idx = NUC.index(ref_poly[i])
        mut_freqs[:, i, :] = np.delete(freq[:, i, :], ref_idx, axis=1)
        mut_counts.append([counts[i][j] for j in range(len(counts[i]))])

    try:
        output_file = os.path.join(out_dir, filename)
        with open(output_file, mode='wb') as f:
            np.savez_compressed(f, mutant_sites=mutant_sites_samp, frequency=mut_freqs, mut_counts=mut_counts, ref_sites=ref_sites)
        print(f"Mutant site frequency file saved to {output_file}")
    except Exception as e:
        print("Error:", e)

def main(args):
    """ determine the maximum frequency that each mutant site reaches in a given region"""
    
    parser = argparse.ArgumentParser(description='Maximum frequency for mutant sites')
    parser.add_argument('-o',             type=str,    default=None,           help='output directory')
    parser.add_argument('--input_file',    type=str,    default=None,           help='input file')
    parser.add_argument('--trim_dir',    type=str,    default=None,           help='un-bp trim dir')
    parser.add_argument('--refFile',      type=str,    default='ref-index.csv', help='reference index file')
    
    arg_list = parser.parse_args(args)
    
    out_dir     = arg_list.o        
    input_file   = arg_list.input_file
    trim_dir = arg_list.trim_dir     
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    print(input_file)
    
    ref_index = pd.read_csv(arg_list.refFile)
    index     = list(ref_index['ref_index'])
    index     = [str(i) for i in index]
    ref_nucs  = list(ref_index['nucleotide'])
    
    find_max_freq(input_file, trim_dir, out_dir, ref_seq=ref_nucs, ref_labels=index)

    #for filename in os.listdir(input_dir):
    #    if filename.endswith(".csv"):
    #        in_file = os.path.join(input_dir, filename)
    #        find_max_freq(in_file, out_dir, ref_seq=ref_nucs, ref_labels=index)

if __name__ == "__main__":
    main(sys.argv[1:])
