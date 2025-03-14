# Libraries, packages

import sys, os

import time

import random

import numpy as np

import pandas as pd

import argparse

def back_proj_format(WRITE_DIR):

    df_dict = {}
    concat_dict = {}

    for root, dirs, files in os.walk(WRITE_DIR):
        subdirectory = os.path.basename(root)
        df_list = []

        for file in files:
            if file.endswith('.csv'):
                path = os.path.join(root, file)
                df_bp = pd.read_csv(path)

                num = int(file.split('[')[1].split(']')[0])
                df_bp['allele'] = num
                df_bp = df_bp[['times', 'rates', 'allele']]

                df_list.append(df_bp)

        df_dict[subdirectory] = df_list

    for subdirectory, df_list in df_dict.items():
        if df_list:
            concat_df = pd.concat(df_list, ignore_index=True)
            concat_dict[subdirectory] = concat_df
            
    return concat_dict

def extract_bp_vec_npz(npz_file, npz_base_name, WRITE_DIR):

    start_tot = time.time()

    concatted = back_proj_format(WRITE_DIR)
    
    bigger_seq = []
    bigger_count = []
    big_times = []
    
    #for subdirectory, concat_df in concatted.items():
    for key, concat_df in concatted.items():
        
        # the unique times for each concatted dataframe corresponding to bp sub-dir
        unique_values_bp = concat_df['times'].unique()
        time_array = [unique_values_bp]
    
        # initialize lists
        big_seq = []
        big_count = []

        # for each time, extract nVec and sVec and populate lists/arrays
        for time_i in unique_values_bp:
            concat_df_vec_pairs = concat_df[concat_df['times']==time_i]
            concat_df_vec_pairs.sort_values(by=['allele'], inplace=True)

            count_array = [x for x in concat_df_vec_pairs['rates']]
            seq_array = [x for x in concat_df_vec_pairs['allele']]

            #seq_array = np.asarray(seq_array)
            seq_array = [[x] for x in seq_array]

            # populate arrays with values
            big_seq.append(seq_array)
            big_count.append(count_array)
            
        # append each simulation value to larger array
        bigger_count.append(big_count)
        bigger_seq.append(big_seq)
        big_times.append(time_array)

    with np.load(npz_file,  allow_pickle=True) as data_mult:
        mutant_sites = data_mult['mutant_sites']
        mutant_sites_all = data_mult['mutant_sites_all']
        n_mutations = data_mult['n_mutations']
        simulations = data_mult['simulations']
        full_nVec = data_mult['full_nVec']
        full_sVec = data_mult['full_sVec']
        times = data_mult['times']
    
    out_str = npz_base_name 

    # save arrays
    f = open(out_str+'-backprojected.npz', mode='w+b')
    np.savez_compressed(f, mutant_sites_all=mutant_sites_all, mutant_sites=mutant_sites, 
                        simulations=simulations, full_nVec=bigger_count, full_sVec=bigger_seq, 
                        n_mutations=n_mutations, times=big_times)
    f.close()

    end_tot = time.time()

    print(f"File saved to {out_str + '-backprojected.npz'} \tTime taken: {(end_tot-start_tot):.03f}s")

    return(f)


def parse_args():
    
    #parse command line arguments passed to script
    parser = argparse.ArgumentParser(
        description='take in back-projected csv files and create npz file formatted for inference')
    
    parser.add_argument('-npz_file', help='file path to stochastic npz file as defined in notebook')
    parser.add_argument('-npz_base_name', help='base name for file naming as defined in notebook')
    parser.add_argument('-WRITE_DIR', help='directory for back-projected estimates')
    return parser.parse_args()

def main():
    #main entry point for module
    args = parse_args()
    npz_file = args.npz_file
    npz_base_name = args.npz_base_name
    WRITE_DIR = args.WRITE_DIR
    
    if not os.path.exists(WRITE_DIR):
        raise Exception('directory does not exist ({0}).'.format(WRITE_DIR))
        
    extract_bp_vec_npz(npz_file, npz_base_name, WRITE_DIR)
    
if __name__ == '__main__':
    main()