import argparse
import os
import pandas as pd
import numpy as np
import sys

def extract_region_name(input_string):
    """
    Extract the region name from the input string.
    """
    # Split the string using '---' as delimiter and take the first part
    region_name = input_string.split('---')[0]
    return region_name

def main(args):
    """ Eliminate sites that don't have more than min_count genomes containing the mutation in the whole time series """
    
    parser = argparse.ArgumentParser(description='Selection coefficients inference')
    parser.add_argument('-o',             type=str,    default=None,           help='output directory')
    parser.add_argument('--mask_dir',      type=str,    default=None,           help='directory containing mask for each site')
    parser.add_argument('--input',        type=str,    default=None,           help='input file containing sequence data')
    parser.add_argument('--refFile',      type=str,    default='ref-index.csv',help='the file containing the reference index and nucleotides')
    parser.add_argument('--gapFile',      type=str,    default='gap-list.npy', help='the file containing sites that have known deletions')
    
    arg_list = parser.parse_args(args)
    
    out_dir   = arg_list.o        
    mask_dir  = arg_list.mask_dir
    file      = arg_list.input
    
    ref_df = pd.read_csv(arg_list.refFile)
    nucs   = np.array(ref_df['nucleotide'])
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
        
    gap_list = np.load(arg_list.gapFile, allow_pickle=True)
    file_tail = os.path.split(file)[-1]
    
    # Extract region name
    region_name = extract_region_name(file_tail)
    mask_path = os.path.join(mask_dir, f'{region_name}-combined-mask.npy')

    print(mask_path)

    main_process(mask_path, file, out_dir, mask_dir, ref_df, nucs, gap_list)

# this is where the mask npz file is accessed
def main_process(mask_path, file, out_dir, mask_dir, ref_df, nucs, gap_list):
    
    file_tail = os.path.split(file)[-1]
    mask_data  = np.load(mask_path, allow_pickle=True)
    seq_data   = get_data(file)
    refs       = np.array(seq_data['ref_sites'])
    muts       = np.array(seq_data['mutant_sites'])
    index  = list(ref_df['ref_index'])
    seq        = seq_data['sequences']
    times      = seq_data['times']
    counts     = seq_data['counts']
    mask       = mask_data

    seqs_new   = np.array([i[mask] for i in seq_data['sequences']])
    muts_new   = muts[mask]
    refs_new   = refs[mask]
    
    idxs = []
    for i in refs_new:
        idx, gap_num = separate_label_idx(str(i))
        idxs.append(int(idx))
    idxs = np.array(idxs)
    mask = np.array([True if (idxs[i] > START_IDX and idxs[i] < END_IDX) else False for i in range(len(idxs))])
    refs_new = refs_new[mask]
    seqs_new = [np.array(i)[mask] for i in seqs_new]

    seqs_new = [''.join(list(i)) for i in seqs_new]
  
    out_file = os.path.join(out_dir, file_tail)
    data = {
        'counts' : seq_data['counts'],
        'times' : times,
        'sequence' : seqs_new
    }
    df = pd.DataFrame(data=data)
    df.to_csv(out_file, index=False) #+ .csv??
        
    sites_data = {'ref_sites' : refs_new}
    df2 = pd.DataFrame(data=sites_data)
    df2.to_csv(out_file[:-4] + '-sites.csv', index=False)

def get_data(file):
    """
    Placeholder function to simulate reading sequence data from a file.
    Replace this with actual data loading logic.
    """
    # Sample data for demonstration purposes
    return {
        'ref_sites': np.array([1, 2, 3]),
        'mutant_sites': np.array([4, 5, 6]),
        'sequences': np.array([['A', 'T', 'G'], ['C', 'G', 'T'], ['G', 'A', 'C']]),
        'times': np.array(['2021-01-01', '2021-02-01', '2021-03-01']),
        'counts': np.array([10, 15, 20])
    }

def separate_label_idx(label):
    """
    Placeholder function to simulate separating label index.
    Replace this with actual label separation logic.
    """
    # Sample logic for demonstration purposes
    idx = label.split('-')[0]
    gap_num = 0
    return idx, gap_num

if __name__ == '__main__': 
    main(sys.argv[1:])
