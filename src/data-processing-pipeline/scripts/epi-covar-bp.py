# coding: utf-8

import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance
import os
import math
import datetime as dt
import subprocess
import pandas as pd
import shutil
import gzip

REF_TAG = 'EPI_ISL_402125'
NUC     = ['-', 'A', 'C', 'G', 'T']

def find_site_index_file(filepath):
    """ Given a sequence file find the correct corresponding file that has the site names"""
    
    directory, file = os.path.split(filepath)
    return filepath[:-4] + '-sites.csv'


def get_data(file, get_seqs=True):
    """ Given a sequence file, get the sequences, dates, and mutant site labels"""
    
    print('sequence file\t', file)
    data = pd.read_csv(file, engine='python', dtype={'sequence': str})
    dates      = list(data['times'])
    counts     = list(data['counts'])
    index_file = find_site_index_file(file)
    index_data = pd.read_csv(index_file)
    ref_sites  = list(index_data['ref_sites'])
    if get_seqs:
        sequences   = np.array([np.array(list(str(i)), dtype=int) for i in list(data['sequence'])])
        seq_lens, seq_counts = np.unique([len(i) for i in sequences], return_counts=True)
        print(f'sequence lengths are {seq_lens}')
        print(f'number of sequences with each length {seq_counts}')
        common_len = seq_lens[np.argmax(seq_counts)]
        mask       = [i for i in range(len(sequences)) if len(sequences[i])==common_len]
        dates      = list(np.array(dates)[mask])
        counts     = list(np.array(counts)[mask])
        sequences  = sequences[mask]
        assert(len(ref_sites)==common_len)
    else:
        sequences = None
    dic = {
        'ref_sites' : ref_sites,
        'sequences' : sequences,
        'dates' : dates,
        'counts' : counts
    }
    return dic


'''def allele_counter(df, d=5):
    """ Calculates the counts for each allele at each time. """
    
    L = len(df['sequences'].iloc[0])
    single = np.zeros((len(np.unique(df['times'])), L * d))
    for df_iter, df_entry in df.iterrows():
        seq = np.array(list(df_entry['sequences']))
        for i in range(len(seq)):
            single[df_iter, (i * d) + seq[i]] += df_entry['counts']
    return single'''



def allele_counter(df, d=5):
    """ Calculates the counts for each allele at each time. """
    
    L = len(df['sequences'].iloc[0])
    unique_times = np.unique(df['times'])
    single = np.zeros((len(unique_times), L * d))
    
    for time_iter, time_value in enumerate(unique_times):
        time_entries = df[df['times'] == time_value]
        
        for _, entry in time_entries.iterrows():
            seq = np.array(list(entry['sequences']))
            for i in range(len(seq)):
                single[time_iter, (i * d) + seq[i]] += entry['counts']
    print(single)
    
    return single


def write_seq_file(seq_file, df):
    f = open(seq_file, 'w')
    print(f)
    for df_iter, df_row in df.iterrows():
        seq    = [str(i) for i in df_row['sequences']]
        counts = df_row['counts']
        #print(counts)
        time   = df_row['times']
        f.write('%d\t%f\t%s\n' % (time, counts, ' '.join(seq)))
    f.close()

def run_mpl(N, q, seq_file, covar_dir, location, counts=False):
    """ Run C++ code that calculates the covariance"""

    covar_file  = f'covar-{location}.dat'
    num_file    = f'num-{location}.dat'
    out_file    = f'out-{location}.dat'
    double_file = f'double-{location}.dat'
    if os.path.exists(os.path.join(covar_dir, covar_file)):
        os.remove(os.path.join(covar_dir, covar_file))
        
    if counts:
        proc = subprocess.Popen(['/Users/liz/Documents/back_projection_home/backprojection-SARS-CoV-2-main/data-processing-pipeline/inference-c++/bin/mpl', '-d', covar_dir, '-i', seq_file, '-o', out_file, '-g', str(0), '-N', f'{N}', '-sc', covar_file, '-sn', num_file, '-q', str(q), '-dc', double_file], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    else: 
        proc = subprocess.Popen(['/Users/liz/Documents/back_projection_home/backprojection-SARS-CoV-2-main/data-processing-pipeline/inference-c++/bin/mpl', '-d', covar_dir, '-i', seq_file, '-o', out_file, '-g', str(0), '-N', f'{N}', '-sc', covar_file, '-sn', num_file, '-q', str(q)], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        print(proc)
    stdout, stderr = proc.communicate()
    return_code    = proc.returncode
    proc.wait()
    if return_code != 0:
        print('directory\t', covar_dir)
        print('file\t', seq_file)
        print(os.listdir(covar_dir))
    return stdout, stderr, return_code


def read_covariance(covar_path, time_varying=False, counts=False, tv_dir=None, location=None, dates_covar=None, coefficient=None):
    """ Read in the covariance file. 
    tv_dir, location, and dates_covar are needed to save the covariance at each time"""

    covar_int     = []
    counter       = 0
    print(covar_path) 
    for line in open(os.path.join(covar_path)):
        if counts:
            covar_int.append(np.array(np.array(line.split(' ')[:-1])[::2], dtype=float))
        else:
            new_line     = line.split(' ')
            new_line[-1] = new_line[-1][:-1]
            covar_int.append(np.array(new_line, dtype=float))
            counter += 1
    return np.array(covar_int, dtype=float)

'''def freqs_from_counts(freq, num_seqs, window=5, min_seqs=200, hard_window=False):
    """ Smooths the counts and then finds the total frequencies"""
    # Determine window size, requiring min_seqs sequences in the beginning and end or a maximum of window days
    max_window  = window
    start_days  = 1
    end_days    = 1
    start_count = num_seqs[0]
    end_count   = num_seqs[-1]
    print(freq)
    while start_count < min_seqs and start_days < window and start_days < len(freq) / 2:
        start_count += num_seqs[start_days]
        start_days  += 1
    while end_count < min_seqs and end_days < window and end_days < len(freq) / 2:
        end_count   += num_seqs[-1-end_days]
        end_days    += 1
    if start_days < max_window and end_days < max_window:
        window = max([start_days, end_days])    
    if window > len(freq / 2):
        window = int(len(freq / 2))
    if window > max_window:
        window = max_window
        
    # if hard_window is True, don't adjust the window size based on the number of sequences in a region
    if hard_window:
        window=max_window
    
    #print(np.shape(freq))
    ret = np.cumsum(freq, axis=0)
    #print(np.shape(ret))
    print('window', window)
    ret[window:] = ret[window:] - ret[:-window]
    result = ret[window - 1:]
        
    # Find total number of sequences in each window from num_seqs
    n_seqs = np.zeros(len(result))
    for i in range(len(result)):
        print(i)
        n_temp = 0
        for j in range(window):
            print(j)
            n_temp += num_seqs[i + j]
        n_seqs[i] += n_temp
            
    #print(len(n_seqs))
    #print(np.shape(result))
    final  = []
    for i in range(len(result)):
        if n_seqs[i]==0:
            final.append(np.zeros(len(result[i])))
        else:
            final.append(result[i] / np.array(n_seqs[i]))
    final  = np.array(final)
    return final'''

'''def freqs_from_counts(counts):
    total = np.sum(counts, axis=1)
    print(total)
    print(counts)
    return counts / (total * 5)'''

def freqs_from_counts(counts, L):
    """Computes the total frequencies without smoothing the counts."""

    final = []
    total = np.sum(counts, axis=1)
    print(total)
    print(counts.shape)
    for i in range(len(counts)):
        if total[i] == 0:
            final.append(np.zeros(len(counts[i])))
        else:
            final.append(counts[i] / np.array(total[i]))
    final = np.swapaxes(np.swapaxes(counts, 0, 1) / [np.sum(counts, axis=1) / L], 0, 1)
    print(f'frequency matrix: {final}')

    
    return final


def main(args):
    """Calculate the covariance matrix and the trajectories for a single region from SARS-CoV-2 data"""
    
    # Read in parameters from command line
    parser = argparse.ArgumentParser(description='Selection coefficients inference')
    parser.add_argument('-o',            type=str,    default='parallel',              help='output directory')
    parser.add_argument('-R',            type=float,  default=2,                       help='the basic reproduction number')
    parser.add_argument('-k',            type=float,  default=0.1,                     help='parameter determining shape of distribution of new infected')
    parser.add_argument('-q',            type=int,    default=5,                       help='number of mutant alleles per site')
    parser.add_argument('--mu',          type=float,  default=0.,                      help='the mutation rate')
    parser.add_argument('--data',        type=str,    default=None,                    help='.npz files containing the counts, sequences, times, and mutant_sites for different locations')
    parser.add_argument('--scratch',     type=str,    default='scratch',               help='scratch directory to write temporary files to')
    parser.add_argument('--record',      type=int,    default=1,                       help='number of generations between samples')
    parser.add_argument('--timed',       type=int,    default=0,                       help='if 0, wont print any time information, if 1 will print some information')
    parser.add_argument('--pop_size',    type=int,    default=10000,                   help='the population size')
    parser.add_argument('--ss_pop_size', type=str,    default=None,                    help='.npy file containing the simulation specific population size')
    parser.add_argument('--ss_k',        type=str,    default=None,                    help='.npy file containing the simulation specific k')
    parser.add_argument('--ss_R',        type=str,    default=None,                    help='.npy file containing the simulation specific R')
    parser.add_argument('--ss_record',   type=str,    default=None,                    help='.npy file containing the simulation specific record')
    parser.add_argument('--outflow',     type=str,    default=None,                    help='.npz file containing the counts and sequences that outflow')
    parser.add_argument('--inflow',      type=str,    default=None,                    help='.npz file containing the counts and sequences that inflow, and the locations')
    parser.add_argument('--window',      type=int,    default=15,                      help='the number of days over which to take the moving average')
    parser.add_argument('--delta_t',     type=int,    default=0,                       help='the amount of dates at the beginning and the end of the time series to use to calculate delta_x')
    parser.add_argument('--mask_site',   type=int,    default=None,                    help='the site to mask (change to WT) when inferring the other coefficients, in order to check effect of linkage (write)')
    parser.add_argument('--minSeqs',     type=int,    default=50,                      help='the minimum number of sequences to use at the end at the beginning to calculate delta_x.')
    parser.add_argument('--mask_group',  type=str,    default=None,                    help='the .npy file containing a list of sites (as nucleotide numbers) to mask for the inference')
    parser.add_argument('--remove_sites',type=str,    default=None,                    help='the sites to eliminate when inferring coefficients')
    parser.add_argument('--decay_rate',  type=float,  default=0,                       help='the exponential decay rate used to correct for infection lasting multiple generations')
    parser.add_argument('--nm_popsize',  type=str,    default=None,                    help='.csv file containing the population sizes used to correct for infection lasting multiple generations')
    parser.add_argument('--delay',       type=int,    default=None,                    help='the delay between the newly infected individuals and the individuals that')
    parser.add_argument('--trajectories',   action='store_true', default=False,  help='whether or not to save the trajectories')
    parser.add_argument('--tv_inference',   action='store_true', default=False,  help='whether or not inference is done at every time')
    parser.add_argument('--find_counts',    action='store_true', default=False,  help='whether or not to calculate the single and double site frequencies')
    
    arg_list = parser.parse_args(args)
    
    out_str       = arg_list.o
    mu            = arg_list.mu
    q             = arg_list.q
    record        = arg_list.record
    window        = arg_list.window
    timed         = arg_list.timed
    input_str     = arg_list.data
    mask_site     = arg_list.mask_site
    decay_rate    = arg_list.decay_rate
    delay         = arg_list.delay
    min_seqs      = arg_list.minSeqs
    scratch_dir   = arg_list.scratch
    
    if arg_list.ss_record:
        record = np.load(arg_list.ss_record)     # The time differences between recorded sequences
    else:
        record = arg_list.record * np.ones(1)
    if arg_list.ss_pop_size:
        pop_size = np.load(arg_list.ss_pop_size) 
    else:
        pop_size = arg_list.pop_size 
    if arg_list.ss_R:
        R = np.load(arg_list.ss_R) # The basic reproductive number                   
    else:
        R = arg_list.R
    if arg_list.ss_k:
        k_ss = np.load(arg_list.ss_k) # The dispersion parameter
    else:
        k_ss = arg_list.k
    if arg_list.delta_t == 0 :
        delta_t = window
    else:
        delta_t = arg_list.delta_t
    if arg_list.remove_sites:
        remove_sites = np.load(arg_list.remove_sites)
    
    # creating directories for c++ files
    if not os.path.exists(out_str):
        os.mkdir(out_str)
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)
    identifier = os.path.split(out_str)[-1]
    covar_dir  = os.path.join(scratch_dir, f'{identifier}-sd-covar-dir')  # the directory that the covariance matrices, will be written to
    if not os.path.exists(covar_dir):
        os.mkdir(covar_dir)     
    print('scratch directory:   ', scratch_dir)
    print('covariance directory:', covar_dir)
    print(f'outfile {out_str}')
    
    # Load the inflowing sequence data if it is known
    if arg_list.inflow:
        inflow_data  = np.load(arg_list.inflow, allow_pickle=True) # sequences that flow migrate into the populations
        in_counts    = inflow_data['counts']        # the number of each sequence at each time
        in_sequences = inflow_data['sequences']     # the sequences at each time
        in_locs      = inflow_data['locations']     # the locations that the sequences flow into
        ### The below only works if a constant population size is used and a single region has inflowing population###
        pop_in       = [np.sum(in_counts[0][i] / np.mean(pop_size)) for i in range(len(in_counts[0]))]  # the total number of inflowing sequences at each time
        #print(pop_in)
            
    if timed > 0:
        t_start = timer()
        print("starting")
    
    # Make status file that records the progress of the script
    filepath    = input_str
    location    = os.path.split(filepath)[-1][:-4]
    status_name = f'covar-test-{location}.csv'
    status_file = open(status_name, 'w')
    status_file.close()
    
    def print2(*args):
        """ Print the status of the processing and save it to a file."""
        stat_file = open(status_name, 'a+')
        line      = [str(i) for i in args]
        string    = '\t'.join(line)
        stat_file.write(string+'\n')
        stat_file.close()
        print(string)
    
    # Load data
    data = get_data(filepath)
    times_temp   = data['dates']
    ref_sites    = np.array(data['ref_sites'])
    counts       = np.array(data['counts'])
    sequences    = [list(i) for i in data['sequences']]
    unique_dates = np.unique(times_temp)
    new_data = {
        'times' : times_temp,
        'counts' : counts,
        'sequences' : sequences
    }
    df = pd.DataFrame(data=new_data)
    print2(ref_sites)
    
    location_data = location.split('-')
    timestamps    = location_data[-6:]
    for i in range(len(timestamps)):
        if len(timestamps[i]) == 1:
            timestamps[i] = '0' + timestamps[i]
                
    dates_full  = np.arange(np.amin(times_temp), np.amax(times_temp) + 1)
    nm_popsize  = pop_size
    L           = len(ref_sites)
    ss_record   = record
    coefficient = pop_size * k_ss * R / (R + k_ss)
    print2('number of sites:', L)
        
    # get data for the specific location
    region = location[:-22]
    if   region[-2:] == '--': region = region[:-2]
    elif region[-1]  ==  '-': region = region[:-1]
        
    # Load the sequences due to travel
    if arg_list.inflow:
        if region in list(in_locs):
            ss_incounts = in_counts[list(in_locs).index(region)]     # Contains numbers of each sequence that migrates to the location for each time
            ss_inseqs   = in_sequences[list(in_locs).index(region)]  # Contains containing the sequences migrating to this location for each time
                
    # Mask all sites in the list of sites given
    if arg_list.mask_group:
        ### NEED TO UPDATE TO DEAL WITH MULTIPLE ALLELES AT EACH SITE
        ref_file   = pd.read_csv('ref-index-short.csv')
        index      = list(ref_file['ref_index'])
        ref_seq    = list(ref_file['nucleotide'])
        ref_poly   = []
        for i in range(len(index)):
            if index[i] in ref_sites:
                ref_poly.append(ref_seq[i])
                
        mutants    = [get_label_new(i) for i in ref_sites]
        mask_group = np.load(arg_list.mask_group, allow_pickle=True)
        mask_sites = [i[:-2] for i in mask_group]
        if np.any([i in mutants for i in mask_sites]):
            mask_nucs  = [NUC.index(i[-1]) for i in mask_group]
            mask_nums  = [ref_sites[mutants.index(i)] for i in mask_sites]
            ref_mask   = [NUC.index(i) for i in ref_poly[mask_nums]]
            mask_idxs  = [mutants.index(i) for i in mask_sites]
            for i in range(len(mask_group)):
                if mask_nums[i] in list(ref_sites):
                    for j in range(len(sVec)):
                        for k in range(len(sVec[i])):
                            if sVec[j][k][mask_idxs[i]] == mask_nucs[i]:
                                sVec[j][k][mask_idxs[i]] = ref_mask[i]
    # mask single site if given                            
    if arg_list.mask_site:
        # mask out 'mask_site' if there is one or multiple
        if mask_site in list(ref_sites):
            site_loc  = list(ref_sites).index(mask_site)
            for i in range(len(sVec)):
                for j in range(len(sVec[i])):
                    sVec[i][j][site_loc] = 0
                        
    # process the sequences that flow to this location from elsewhere, if they are known.
    ### DOES THIS NEED TO BE FIXED FOR MULTIPLE STATES AT EACH SITE?
    if arg_list.inflow and not mask_site:
        if region in list(in_locs):
            traj_temp = trajectory_calc(nVec, sVec, ref_sites, d=q)
            print2('number of time points in the region  \t', len(nVec))
            print2('number of time points inflowing      \t', len(ss_incounts))
            print2('number of time points in trajectories\t', len(traj_temp))
            inflow_term = allele_counter_in(ss_inseqs, ss_incounts, ref_sites, traj_temp, pop_in, pop_size, k_ss, R, len(nVec)) # Calculates the correction due to migration
            print2(f'region {region} present')
        else:
            inflow_term = np.zeros(len(ref_sites) * q)
    else:
        inflow_term = np.zeros(len(ref_sites) * q)
    
    ### Run inference trimming the data by submission date
    if timed > 0:
        bootstrap_start_time = timer()
            
    df = df.sort_values(by='times')
    print(df)
        
    # Write nVec and sVec to file
    seq_file = f'seqs-{location}.dat'
    seq_path = os.path.join(covar_dir, seq_file)
    write_seq_file(seq_path, df)
    
    if arg_list.find_counts:
        stdout, stderr, exit_code = run_mpl(pop_size, q, seq_file, covar_dir, location, counts=True)
    else:
        stdout, stderr, exit_code = run_mpl(pop_size, q, seq_file, covar_dir, location, counts=False)
    
    if exit_code != 0:
        print2(exit_code)
        print2(stdout)
        print2(stderr)
        print2("mistake")
    
    covar_file = f'covar-{location}.dat'
    if arg_list.tv_inference:
        cov_path_old = os.path.join(covar_dir, covar_file)
        proc         = subprocess.Popen(['gzip', cov_path_old])
        proc.wait()
        if os.path.exists(cov_path_old):
            os.remove(cov_path_old)
        
    # read in covariance file  
    if arg_list.tv_inference:
        tv_covar_dir     = os.path.join(out_str, 'tv_covar')
        t_unique, counts = np.unique(times_temp, return_counts=True)
        total_seqs    = []
        for i in range(len(dates_full)):
            if dates_full[i] in t_unique:
                total_seqs.append(counts[list(t_unique).index(dates_full[i])])
            else: 
                total_seqs.append(0)
        covar_path = os.path.join(covar_dir, covar_file + '.gz')
        dates_covar   = [dates_full[i] for i in range(1, len(dates_full)) if total_seqs[i-1]!=0]
        dates_nocovar = [dates_full[i] for i in range(1, len(dates_full)) if (dates_full[i] not in dates_covar and dates_full[i]>dates_covar[0])]
        covar_int     = read_covariance(covar_path, time_varying=True, tv_dir=tv_covar_dir, 
                                        location=location, dates_covar=dates_covar, coefficient=coefficient)
        print2(f'dates without a covariance matrix: {dates_nocovar}')
        for t in dates_nocovar:
            old_file = os.path.join(tv_covar_dir, location + f'___{t-1}.npz')
            new_file = os.path.join(tv_covar_dir, location + f'___{t}.npz')
            shutil.copy(old_file, new_file)
        covar_int = []
    elif arg_list.find_counts:
        double_file = f'double-{location}.dat'
        covar_path  = os.path.join(covar_dir, double_file)
        covar_int   = read_covariance(covar_path, counts=True)
    else:
        covar_path = os.path.join(covar_dir, covar_file)
        print(covar_path)
        covar_int  = read_covariance(covar_path) * coefficient / 5 
    
    # delete files
    for file in [seq_file, f'covar-{location}.dat', f'num-{location}.dat', f'out-{location}.dat', f'covar-{location}.dat.gz']:
        if os.path.exists(os.path.join(covar_dir, file)):
            os.remove(os.path.join(covar_dir, file))
    print2(ref_sites)
    
    # calculate the frequencies
    counts = allele_counter(df, d=q)
    print(f'counts matrix: {counts}')
    unique_dates = np.unique(list(df['times']))
    num_seqs     = [np.sum(df[df['times']==i]['counts']) for i in unique_dates]
    #total_seqs    = []
    #for i in range(len(dates_full)):
    #    if dates_full[i] in unique_dates:
    #        total_seqs.append(num_seqs[list(unique_dates).index(dates_full[i])])
            #print(total_seqs)
    #    else: 
    #        total_seqs.append(0)
    #if arg_list.tv_inference and 'united kingdom' not in location:
    #    counts = freqs_from_counts(counts, total_seqs, min_seqs=min_seqs)
    if not arg_list.find_counts:
        counts = freqs_from_counts(counts, L)
    else:
        counts = np.sum(counts, axis=0)
        
    # save data
    file = os.path.join(out_str, location + '.npz')
    g = open(file, mode='wb')
    if not arg_list.find_counts:
        np.savez_compressed(
            g, 
            location=location, 
            times=dates_full[-len(counts):],
            times_full=dates_full,
            ref_sites=ref_sites, 
            allele_number=ref_sites,  
            k=k_ss, N=pop_size, R=R,
            covar=covar_int, 
            counts=counts, 
            inflow=inflow_term
        )
    else:
        np.savez_compressed(
            g, 
            location=location, 
            ref_sites=ref_sites, 
            allele_number=ref_sites,  
            double_counts=covar_int, 
            single_counts=counts, 
        )
    g.close()                    
        
    if timed > 0:
        new_run_time = timer()
        print2(f'{new_run_time - bootstrap_start_time} seconds to run the covariance')
        bootstrap_start_time = new_run_time

    os.rmdir(covar_dir)   


if __name__ == '__main__': 
    main(sys.argv[1:])

