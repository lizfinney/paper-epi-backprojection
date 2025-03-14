import sys
import argparse
import os
import csv
import numpy as np
import pandas as pd
from itertools import islice
import datetime



def in_formatted(in_dir, file, thresh=0.04):

    #df = pd.read_csv("africa-kenya-formatted.csv")
    df = pd.read_csv(file)
    df = pd.read_csv(os.path.join(in_dir, file))

    # exclude the sequence column in the first column
    for col in df.columns[1:]:
        # get rid of values less than 4% (based on prob distribution value)
        # basically truncating the probability distribution
        df[col] = df[col].apply(lambda x: 0 if pd.to_numeric(x, errors='coerce') < thresh else x)

    # get rid of columns that now have all 0s after counts are scrubbed below an x% threshold (default 4%)
    df = df.loc[:, (df != 0).any(axis=0)]

    # list of columns, pop off sequence column name
    dflist=list(df)
    dflist.pop(0)

    # integer values
    date_list = [int(i) for i in dflist]

    return df, date_list


# function finds windows of consecutive time values, default n=5 days for window

def consecutive_windows(seq, n=5):
    it = iter(seq)
    result = tuple(islice(it, n))
    
    while len(result) == n:
        if all(result[i] == result[i-1] + 1 for i in range(1, n)):
            yield result
        try:
            result = result[1:] + (next(it),)
        except StopIteration:
            break
            
# write below as a function that takes in different values for filtering criteria


def threshold_counts(in_dir, file, x=10, thresh=0.04, n=5):

    window_collection = []

    ''' after consecutive days found, checks consecutive days to see if sum of these date columns is a value above the desired threshold 
        paper did 20 sequences in 5 day window
        lowered this threshold to capture the timespans that originally plummeted to 0 and then went back up, so 10 is default
        update: keep 20 at default to avoid memory issues downstream
    '''
    
    df, date_list = in_formatted(in_dir, file, thresh)

    for window_result in consecutive_windows(date_list, n):
        total_window = []
        for i in window_result:
            total = df[f'{i}'].sum()
            total_window.append(total)
        total_sum = sum(total_window)
        if total_sum > x:
            window_collection.append(window_result)

    # extract unique integer values from window lists
    unique_values = set()
    for window in window_collection:
        unique_values.update(window)

    unique_values = sorted(unique_values)

    return(unique_values)


def check_consecutive_groups(data):
    consecutive_groups = []
    current_group = []
    labeled_groups = {}
    group_counter = 1


    # for each unique date value in the windows that pass count threshold, check to see
    # if it is part of a series of consecutive dates. Previous code did 20 day trimming.

    for num in data:
        if not current_group or num == current_group[-1] + 1:
            current_group.append(num)
        else:
            if len(current_group) >= 20:
                label = f"group_{group_counter}"
                labeled_groups[label] = [str(i) for i in current_group]
                group_counter += 1
            current_group = [num]

    # check the last group
    if len(current_group) >= 20:
        label = f"group_{group_counter}"
        labeled_groups[label] = [str(i) for i in current_group]

    return labeled_groups

def process_trim_intervals(labeled_dict, location, in_dir, out_dir, file, thresh=0.04):

    df, _ = in_formatted(in_dir, file, thresh)
    #print(df)

    for label, group in labeled_dict.items():
        print(label)
        print(group)
        

        keep = labeled_dict.get(label)
        print(keep)
        #keep = labeled_dict.get(f'group_{group}')

        # subsetted dataframes

        df_subset = df[keep]

        # add sequence column back on

        df_main = pd.concat([df['sequence'], df_subset], axis=1)
 
        # identify rows where all values after the sequence column are all 0
        # (sequence column is never 0)
        df_filter = df_main[df_main.iloc[:, 1:].ne(0).any(axis=1)]

        ''' this part puts the data back into a form that is closer to what the script needs to run
            where sequences are listed multiple times for each date instance.
            This is not ideal but until there are memory issues this is the easiest thing to do.
        '''

        # resets the below index so that this actually works (do I need to keep the index?)

        df_reset = df_filter.reset_index(drop=True)

        # list objects for new dataframe

        seq = []
        times = []
        counts = []

        # iterate through the rows and check non-zero values

        for index, row in df_reset.iterrows():
            non_zero_values = []
            first_col_value = row['sequence']
            for col, val in row.items():
                
                # if column is not the sequence column and the value in the column is not 0
                if col != 'sequence' and val != 0:
                    non_zero_values.append((val, col))
            
            # i looks like (count, time) multiplied for how many instances there are for each sequence
            for i in non_zero_values:
                # add sequence for each date/count object
                seq.append(first_col_value)
                counts.append(i[0])
                times.append(int(i[1]))
                                 #first_col_value, i)
            #print(f"For row {index}, non-zero values: {non_zero_values}")

        df_trim = pd.DataFrame({
            'sequence': seq,
            'times': times,
            'counts': counts
        })

        trimtimes = df_trim['times'].unique()

        # file naming in terms of date ranges

        location = str(location)
        ref_date = datetime.date(2020, 1, 1)

        start_year, end_year   = (ref_date + datetime.timedelta(int(trimtimes[0]))).year,   (ref_date + datetime.timedelta(int(trimtimes[-1]))).year
        start_month, end_month = (ref_date + datetime.timedelta(int(trimtimes[0]))).month,  (ref_date + datetime.timedelta(int(trimtimes[-1]))).month
        start_day, end_day     = (ref_date + datetime.timedelta(int(trimtimes[0]))).day,    (ref_date + datetime.timedelta(int(trimtimes[-1]))).day
        out_file = os.path.join(out_dir, f'{location}---{start_year}-{start_month}-{start_day}-{end_year}-{end_month}-{end_day}.csv')

        print(f'out file is {out_file}')

        df_trim.to_csv(out_file, index=False)


def main(args):
    """ Trims times in the regional sequence data after deconvolution method is completed and formatted """
    
    parser = argparse.ArgumentParser(description='deconvolved regional time trimming')
    parser.add_argument('--input_dir',     type=str,    default=None,     help='input directory')
    parser.add_argument('-o',              type=str,    default=None,     help='output directory')
    parser.add_argument('--location',       type=str,    default=None,     help='region')
    parser.add_argument('-x',              type=int,    default=20,       help='minimum counts')
    parser.add_argument('--thresh',         type=float,  default=0.04,     help='minimum single count value')
    parser.add_argument('-n',              type=int,    default=5,        help='window length')
    
    arg_list = parser.parse_args(args)
            
    input_dir   = arg_list.input_dir
    out_dir  = arg_list.o
    location = arg_list.location
    x        = arg_list.x
    thresh   = arg_list.thresh
    n        = arg_list.n

    file = os.path.join(input_dir, f'{location}-formatted.csv')
    
    print(input_dir)
    print(file)

    unique_values = threshold_counts(input_dir, file, x, thresh, n)
    result = check_consecutive_groups(unique_values)
    process_trim_intervals(result, location, input_dir, out_dir, file, thresh)
        
    print(f'saved to {out_dir}')

    
if __name__ == '__main__': 
    main(sys.argv[1:])
