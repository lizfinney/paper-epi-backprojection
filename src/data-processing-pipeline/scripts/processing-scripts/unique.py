import os
import pandas as pd
import argparse
import sys

def process_unq(input_csv, scratch_out, export_out, root_dir):

    seq_df = pd.read_csv(input_csv)

    # throw out 0th row which is the reference sequence
    #seq_df.drop(index=seq_df.index[0], axis=0, inplace=True)

    # factorizes, creating new seq_code column of factorization values

    seq_df['seq_code'] = pd.factorize(seq_df['sequence'])[0]


    # creates subset dataframe of duplicated values, adds to list

    df_sub = seq_df.loc[seq_df['seq_code'].duplicated(), :]
    a_list = df_sub['seq_code'].tolist()

    # creates unique sequence dataframe and nonunique sequence dataframe

    seq_df_unique = seq_df[~seq_df['seq_code'].isin(a_list)]
    seq_df_nonunique = seq_df[seq_df['seq_code'].isin(a_list)]

    # throw seq code out, reindex, loop through

    nonunique_subset = seq_df_nonunique[["accession", "date", "submission_date", "sequence"]]
    nonunique_subset = nonunique_subset.reset_index(drop=True)

    # refactors seq_code column

    nonunique_subset['seq_code'] = pd.factorize(nonunique_subset['sequence'])[0]

    # groups with extra 'junk' columns, possibly needed to put everything back together

    grouped_scratch = nonunique_subset.groupby(nonunique_subset.seq_code)

    # write groups to separate csv files

    for seq_code, data in grouped_scratch:
        data.to_csv(f"{scratch_out}/{seq_code}.csv")
    
    # groups without extra junk columns
    
    date_count = nonunique_subset[["date","seq_code"]]
    grouped_out = date_count.groupby(date_count.seq_code)
    
    for seq_code, data in grouped_out:
        data_counts = pd.DataFrame(data['date'].value_counts().reset_index())
        data_counts.columns = ['times', 'counts']
        data_counts.to_csv(f"{export_out}/{seq_code}.csv", index=False)


        #data3 = data3.rename(columns = {'index':'times', 'date':'counts'})

    # make unique sequence file
    
    unique_subset = seq_df_unique[["date", "sequence"]]
    unique_subset = unique_subset.rename(columns = {'date':'times', 'sequence':'sequences'})
    unique_subset.to_csv(os.path.join(root_dir, 'unique-sequences.csv'), index=False)
    #unique_subset.to_csv(f"{export_out}/unique-sequences.csv", index = False)

    # creating index file
    
    l = []
    s = []

    for f in os.listdir(scratch_out):
        g = f.strip('.csv')
        l.append(int(g))
        seq_file = pd.read_csv(os.path.join(scratch_out, f))
        seq_df = pd.DataFrame(seq_file)
        #print(seq_df)
        s.append(seq_df['sequence'][0])
    

    index_seq_list = pd.DataFrame(
        {'index': l,
         'sequences': s
        })
    index_seq_list
    new_index = index_seq_list.sort_values(
                                    by='index',
                                    ascending=True)
    
    final = new_index.reset_index()[['index', 'sequences']]
    final.to_csv(os.path.join(root_dir, 'sequence-index.csv'), index=False)
    #final.to_csv(f"{export_out}/sequence-index.csv", index=False)

def main(args):
    ''' formats csv into output folders read by back-projection for nonunique seqs, creates file for unique-seqs'''

    parser = argparse.ArgumentParser(description= 'processing into unique sequences and nonunique sequences')
    parser.add_argument('--input',         type=str,     default=None,     help='input sequence csv file')
    parser.add_argument('-o',              type=str,     default=None,     help='out directory for bp formatted folders')

    arg_list = parser.parse_args(args)

    input_file         = arg_list.input
    out                = arg_list.o

    file_tail = os.path.split(input_file)[-1]
    file_root = os.path.splitext(file_tail)[0]

    root_dir = os.path.join(out, file_root)
    if not os.path.exists(root_dir):
        os.mkdir(root_dir)

    sub_scratch = os.path.join(root_dir, f'scratch-{file_root}')
    sub_out = os.path.join(root_dir, f'output-{file_root}')

    if not os.path.exists(sub_scratch):
        os.mkdir(sub_scratch)

    if not os.path.exists(sub_out):
        os.mkdir(sub_out)

    
    process_unq(input_file, sub_scratch, sub_out, root_dir)

if __name__ == "__main__":
    main(sys.argv[1:])
