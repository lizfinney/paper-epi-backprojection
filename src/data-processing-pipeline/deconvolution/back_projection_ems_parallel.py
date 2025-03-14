# Libraries, packages

from importlib import reload
import sys, os
from copy import deepcopy

import numpy as np

import scipy as sp
import scipy.stats as st

import pandas as pd

import matplotlib
from matplotlib import cm
import matplotlib.pyplot as plt

import argparse

# Global Variables

# EM algorithm

def get_em_lambda(f,   # incubation time distribution, a data frame with columns dt and weight
                  Y,   # observations, a data frame with columns times and counts
                  lam  # infection rates, a data frame with columns times and rates
                  ):
    
    lam_star = lam.copy(deep=True)
    times = lam['times']
    times_dt = f['dt']
    
    for idx in range(len(lam)):
        
        # time for the current lambda value
        t_lam = lam.loc[idx, 'times']
        
        lam_temp = 0
        f_norm   = 0
        for dt in times_dt:
            if (t_lam+dt) not in Y['times'].values:
                continue
            
            Y_tdt = Y[Y['times']==t_lam+dt]['counts'].values[0]  # observations at t + dt
            lam_t = lam[lam['times']==t_lam]['rates'].values[0]  # rate at t
            f_dt  = f[f['dt']==dt]['weight'].values[0]           # incubation prob at dt
            
            #print(idx, t_lam, dt, Y_tdt, lam_t, f_dt)
            
            f_norm += f_dt
            
            # no contribution to lambda* when any of above are zero
            num = Y_tdt * lam_t * f_dt
            if num==0:
                continue
            
            denom = 0
            for t in times:
                ddt = t_lam+dt-t
                if ddt in f['dt']:
                    denom += lam[lam['times']==t]['rates'].values[0] * f[f['dt']==ddt]['weight'].values[0]
    
            lam_temp += num/denom
        
        if f_norm>0:
            lam_star.loc[idx, 'rates'] = lam_temp/f_norm
            
        else:
            lam_star.loc[idx, 'rates'] = 0
        
    return lam_star
                                                                         
# smoothing, S step

def get_smooth_lambda(w,   # smoothing function, a data frame with columns dt and weight
                      lam  # infection rates, a data frame with columns times and rates
                      ):
    
    lam_smooth = lam.copy(deep=True)
    times = lam['times']
    times_dt = w['dt']
    
    for idx in range(len(lam)):
    
        # time for the current lambda value
        t_lam = lam.loc[idx, 'times']
        
        lam_temp = 0
        w_norm   = 0
        for dt in times_dt:
            if (t_lam+dt) not in lam['times'].values:
                continue
            
            lam_tdt = lam[lam['times']==t_lam+dt]['rates'].values[0]  # rate at t + dt
            w_dt    = w[w['dt']==dt]['weight'].values[0]              # smoothing weight at dt
            
            lam_temp += lam_tdt * w_dt
            w_norm   += w_dt
            
        if w_norm>0:
            lam_smooth.loc[idx, 'rates'] = lam_temp/w_norm
        
        else: 
            lam_smooth.loc[idx, 'rates'] = 0
        
    return lam_smooth
                                                                         
                                                                         
#Reading in data and quality checks (no array overruns, etc.)
def back_projection_fit_parallel(input_file, write_dir, incubation_file, smoothing_file):
    
    #incubation_file = incubation_file
    #smoothing_file = smoothing_file
    #observation_files = [os.path.join(read_dir, f) for f in os.listdir(read_dir) if f.endswith(".csv")]
            
    #for file in observation_files:

    df_inc = pd.read_csv(incubation_file, memory_map=True)
    df_smooth = pd.read_csv(smoothing_file, memory_map=True)
    df_obs = pd.read_csv(file, memory_map=True)

    # expand indices such that observations (and lambda extend +/- dt around current boundaries)
    dt_min = np.min([np.min(df_inc['dt']), 0]) + np.min([np.min(df_smooth['dt']), 0])
    dt_max = np.max([np.max(df_inc['dt']), 0]) + np.max([np.max(df_smooth['dt']), 0])

    t_min = np.min(df_obs['times']) + dt_min - 15
    t_max = np.max(df_obs['times']) + dt_max + 1

    for t in range(t_min, t_max+1):
        if t not in df_obs['times'].values:
            row = dict(times=t, counts=0)
            df_obs = pd.concat([df_obs, pd.DataFrame([row])], ignore_index=True)


    # copy lambda data frame, uniform initialization
    df_lam = df_obs.copy(deep=True)
    df_lam['rates'] = 1/float(len(df_lam))


    # sanity checking
    if np.min(df_smooth['weight'])<0:
        print('Smoothing function contains negative values!')

    if np.min(df_inc['weight'])<0:
        print('Incubation time distribution contains negative values!')

    if np.min(df_obs['counts'])<0:
        print('Counts file contains negative values!')

    times = df_lam['times'].values
    for t in times:
        n_rows = len(df_lam[df_lam['times']==t])
        if n_rows!=1:
            print('Got %d observation rows with time=%d!' % (n_rows, t))

    # convergence conditions
    max_iter = 250
    min_dlam = 1e-4

    # initial conditions
    df_lam_init = df_lam.copy(deep=True)
    df_lam_last = df_lam.copy(deep=True)
    curr_iter   = 0

    print('iter\tepsilon')

    while True:
        
        curr_iter += 1
        
        df_lam = get_em_lambda(f=df_inc, Y=df_obs, lam=df_lam_last)
        df_lam = get_smooth_lambda(w=df_smooth, lam=df_lam)

        curr_lam = np.array(df_lam['rates'])
        last_lam = np.array(df_lam_last['rates'])

        eps = np.sum(np.absolute(curr_lam - last_lam)) / np.sum(last_lam)

        print('%d\t%.2e' % (curr_iter, eps))

        df_lam_last = df_lam

        if (curr_iter>=max_iter) or (eps<=min_dlam):
            break

    df_lam_last.sort_values(by=['times'], inplace=True)
    df_lam_last = df_lam_last[df_lam_last['rates'] != 0]

    output_file = os.path.join(write_dir, os.path.basename(file.split('.')[0] + '_rates.csv'))

    #output_base = os.path.splitext(output_file)[0]  # Get the base name of the output file
    #output_rates = output_base + '_rates.csv'  # Append "_rates" to the base name
    # Write result to output file
    #df_lam_last.to_csv(output_rates, index=False)
    df_lam_last.to_csv(output_file)
                                                                         
def parse_args():
    #parse command line arguments passed to script
    parser = argparse.ArgumentParser(
        description='run backprojection and write output')
    
    parser.add_argument('-input_file', help='read csv')
    parser.add_argument('-write_dir', help='write directory')
    parser.add_argument('-incubation_file', help='incubation distribution csv file')
    parser.add_argument('-smoothing_file', help='binomial smoothing weights csv file')
    
    
    return parser.parse_args()

def main():
    #main entry point for module
    args = parse_args()
    input_file = args.input_file
    write_dir = args.write_dir
    incubation_file = args.incubation_file
    smoothing_file = args.smoothing_file
    
    if not os.path.exists(input_file):
        raise FileNotFoundError('Input file does not exist: {}'.format(input_file))

        
    back_projection_fit_parallel(input_file, write_dir, incubation_file, smoothing_file)
    
if __name__ == '__main__':
    main()