import pandas as pd
import datetime as dt 
import argparse
from tqdm import tqdm

METADATA_COLS    = [   'accession', 'virus_name',            'date', 'location',             'location_additional', 'submission_date']
METADATA_OCOLS   = ['Accession ID', 'Virus name', 'Collection date', 'Location', 'Additional location information', 'Submission date']
METADATA_XFORM = [lambda x: str(x), lambda x: str(x), lambda x: str(x), 
                  lambda x: str(x), lambda x: str(x), lambda x: str(x), 
                  lambda x: str(x)]

def process_tsv(input_file, output_file):
    accessions = set()
    header = ','.join(METADATA_COLS)
    use_cols = ['Accession ID', 'Virus name', 'Collection date', 'Location', 
            'Additional location information', 'Submission date', 'Host']
    
    df = pd.read_csv(input_file, sep='\t', usecols=use_cols)
    f = open(output_file, 'w')
    f.write('%s\n' % header)
    
    for df_iter, df_entry in tqdm(df.iterrows(), total=df.shape[0]):
        entry = []
        unique_id = False
        valid_date = False
        
        # check that host is human
        if 'Host' in df_entry:
            if df_entry.Host != 'Human' and df_entry.Host != 'human':
                continue
        
        # get accession
        acc_idx = METADATA_COLS.index('accession')
        if 'Accession ID' not in df_entry:
            continue
        xformed = METADATA_XFORM[acc_idx](df_entry[METADATA_OCOLS[acc_idx]])
        acc = METADATA_XFORM[acc_idx](df_entry[METADATA_OCOLS[acc_idx]])

        if xformed not in accessions:
            accessions.add(xformed)
            entry.append(xformed)
            unique_id = True
            
        # get virus name 
        if 'Virus name' in df_entry:
            name_idx = METADATA_COLS.index('virus_name')
            xformed = METADATA_XFORM[name_idx](df_entry[METADATA_OCOLS[name_idx]])
            xformed = xformed.replace(',', ' ')
            entry.append(xformed)
        else:
            entry.append('NA')
            
        # get date
        if unique_id:
            date_idx = METADATA_COLS.index('date')
            xformed = METADATA_XFORM[date_idx](df_entry[METADATA_OCOLS[date_idx]])

            try:
                dt.date.fromisoformat(xformed)
                valid_date = True
            except:
                valid_date = False
                continue
                
            entry.append(xformed)
            
        # get location 
        if unique_id and valid_date:
            loc_idx = METADATA_COLS.index('location')
            xformed = METADATA_XFORM[loc_idx](df_entry[METADATA_OCOLS[loc_idx]])
            xformed = xformed.replace(',', ' ')
            entry.append(xformed)
        
            loc_add_idx = METADATA_COLS.index('location_additional')
            xformed = METADATA_XFORM[loc_add_idx](str(df_entry[METADATA_OCOLS[loc_add_idx]]))
            xformed = xformed.replace(',', ' ')
            entry.append(xformed)
        
        # get submission date
        if unique_id and valid_date:
            sub_idx = METADATA_COLS.index('submission_date')
            xformed = METADATA_XFORM[sub_idx](df_entry[METADATA_OCOLS[sub_idx]])

            try:
                dt.date.fromisoformat(xformed)
                valid_sub_date = True
            except:
                valid_sub_date = False
                continue
            
            entry.append(xformed)
        
        # save data
        if unique_id and valid_date and valid_sub_date:
            f.write('%s\n' % ','.join(entry))
            
    # Add reference if not in the list of accessions
    REF_TAG = 'EPI_ISL_402125'
    if REF_TAG not in accessions:
        entry = []
        for i in range(len(METADATA_COLS)):
            if METADATA_COLS[i] == 'Accession ID':
                entry.append(REF_TAG)
            elif METADATA_COLS[i] == 'Collection date' or METADATA_COLS[i] == 'Submission date':
                entry.append('2020-01-01')
            else:
                entry.append('reference')
        f.write('%s\n' % ','.join(entry))
        
    f.close()
    
def main():
    parser = argparse.ArgumentParser(description='Process TSV file.')
    parser.add_argument('input_file', type=str, help='Path to the input TSV file')
    parser.add_argument('output_file', type=str, help='Path to the output TSV file')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    process_tsv(input_file, output_file)


if __name__ == '__main__':
    main()

