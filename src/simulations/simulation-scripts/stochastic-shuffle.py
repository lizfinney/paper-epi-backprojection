# Libraries, packages

import sys, os

import time

import random

import numpy as np

import pandas as pd

import argparse


# takes in Times, nVec, sVec
# for each sequence at each time, iterates through counts to reassign alleles to scrambled times


def shuffle_times_mult(npz_file):
    
    # load data for checks

    with np.load(npz_file,  allow_pickle=True) as data_mult:
        mutant_sites_all = data_mult['mutant_sites_all']
        simulations = data_mult['simulations']
        count_vec = data_mult['full_nVec']
        seq_vec = data_mult['full_sVec']
        times_vec = data_mult['times']
    
    
    
    shape_dist = 5.807
    scale_dist = 0.948
    # empty dictionary for dataframes
    df_dict = {}

    for i in range(len(count_vec)):

        start = time.time()
        # i is simulation number, if 10 sims, 0 through 9
        # initialize new vectors for each simulation
        times_vec_new = []
        count_vec_new = []
        seq_vec_new = []
        
        for j in range(len(count_vec[i])):
            # j is each time point, from 0 to 100
            for k in range(len(count_vec[i][j])):
                # k is index for each single vector
                count = count_vec[i][j][k]
                allele = seq_vec[i][j][k]
            
                for l in range(count):
                    # index of counts
                    time_l = times_vec[i][j]
                    dist_samp = np.random.gamma(shape=shape_dist, scale=scale_dist, size=1)
                    mean_shifted_value = dist_samp[0] - shape_dist*scale_dist
                    time_shuff = time_l + int(round(mean_shifted_value))
                    #time_shuff = time + int(round(dist_samp[0]))
                    times_vec_new.append(time_shuff)
                    count_vec_new.append(l)
                    seq_vec_new.append(allele)
    
        # dictionary of lists, write to dataframe
        vec_table = {'times': times_vec_new, 'count_index': count_vec_new, 'allele': seq_vec_new}
        df_dict[f'df_{i}'] = pd.DataFrame(vec_table)
        
        new_dict = {}
        
        for key, df in df_dict.items():
            df['allele'] = df['allele'].astype(str)
        
            unique_pairs = df.groupby(['times', 'allele']).size().reset_index(name='counts')
            new_dict[key] = unique_pairs

        end = time.time()

        print(f"Simulation: {i}\tTime taken: {(end-start)*10**3:.03f}ms")
        
    return new_dict


def extract_shuffled_vec_npz(npz_file, npz_base_name, PROJ_DIR):

    start_tot = time.time()

    # Create the output folder if it doesn't exist
    os.makedirs(PROJ_DIR, exist_ok=True)
    
    bigger_seq = []
    bigger_count = []
    big_times = []
    
    df_dict = shuffle_times_mult(npz_file)

    for key, df in df_dict.items():

        # create sub directories for back-projection
        SUB_DIR = os.path.join(PROJ_DIR, f'{key}')
        os.makedirs(SUB_DIR, exist_ok=True)

        # the unique times, and unique sequence allele values for each dataframe

        unique_seq = df['allele'].unique()
        unique_values = df['times'].unique()
        time_array = [unique_values]

        # initialize lists
        big_seq = []
        big_count = []

        # for each time, extract nVec and sVec and populate lists/arrays
        for time_i in unique_values:
            df_vec_pairs = df[df['times']==time_i]

            count_array = [x for x in df_vec_pairs['counts']]
            seq_array = [x for x in df_vec_pairs['allele']]

            # fixes formatting issues, stupid but as of yet unavoidable
            seq_array = [i.replace(']' , '') for i in seq_array]
            seq_array = [i.replace('[' , '') for i in seq_array]
            seq_array = [int(x) for x in seq_array]

            seq_array = np.asarray(seq_array)
            seq_array = [[x] for x in seq_array]

            # populate arrays with values
            big_seq.append(seq_array)
            big_count.append(count_array)
            

        # append each simulation value to larger array
        bigger_count.append(big_count)
        bigger_seq.append(big_seq)
        big_times.append(time_array)

        for allele in unique_seq:
            df_allele = df[df['allele'] == allele]
            df_allele = df_allele.drop(columns=['allele'])

            # Define the output filename based on the allele value
            output_filename = os.path.join(SUB_DIR, f'dataframe_{allele}.csv')
            df_allele.to_csv(output_filename, index=False)

            print(f"Saved dataframe_{allele} to the {key} sub-directory")

        
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
    f = open(out_str+'-shuffled.npz', mode='w+b')
    np.savez_compressed(f, mutant_sites_all=mutant_sites_all, mutant_sites=mutant_sites, 
                        simulations=simulations, full_nVec=bigger_count, full_sVec=bigger_seq, 
                        n_mutations=n_mutations, times=big_times)
    f.close()

    end_tot = time.time()

    print(f"File saved to {out_str + '-shuffled.npz'} \tTime taken: {(end_tot-start_tot):.03f}s")
    
    return(f)


def parse_args():
    
    #parse command line arguments passed to script
    parser = argparse.ArgumentParser(
        description='import stochastic simulation and write time-shuffled output')
    
    parser.add_argument('-npz_file', help='file path to stochastic npz file as defined in notebook')
    parser.add_argument('-npz_base_name', help='base name for file naming as defined in notebook')
    parser.add_argument('-PROJ_DIR', help='directory for back-projecting')
    return parser.parse_args()

def main():
    #main entry point for module
    args = parse_args()
    npz_file = args.npz_file
    npz_base_name = args.npz_base_name
    PROJ_DIR = args.PROJ_DIR
    
    if not os.path.exists(PROJ_DIR):
        raise Exception('directory does not exist ({0}).'.format(PROJ_DIR))
        
    extract_shuffled_vec_npz(npz_file, npz_base_name, PROJ_DIR)
    
if __name__ == '__main__':
    main()