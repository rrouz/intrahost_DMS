#!/usr/bin/env python3

"""Read in the DMS data from the Flu_H5_DMS_data repository and calculate the DMS scores for each sequence"""

import pandas as pd
import argparse
from tqdm import tqdm
import re

def read_dms_data(dms_df, h5_site_header):
    # convert the DMS data to a json object dms_data['positions']['mutation'] = score
    # Get all columns except known metadata columns
    metadata_columns = {'site', 'wildtype', 'mutant', 'antibody_set'}
    columns_to_retrieve = [col for col in dms_df.columns if col not in metadata_columns]
    dms_data = {}
    
    # First handle antibody_set separately
    if 'antibody_set' in dms_df.columns:
        dms_data['antibody_set'] = {}
        for index, row in dms_df.iterrows():
            site = row[h5_site_header]
            mutant = row['mutant']
            antibody_set = row['antibody_set']
            if site not in dms_data['antibody_set']:
                dms_data['antibody_set'][site] = {}
            dms_data['antibody_set'][site][mutant] = antibody_set
    
    # Then handle all other columns
    for column in tqdm(columns_to_retrieve):
        if column != 'antibody_set':  # Skip antibody_set as it's handled separately
            dms_data[column] = {}
            for index, row in dms_df.iterrows():
                site = row[h5_site_header]
                mutant = row['mutant']
                score = row[column]
                if site not in dms_data[column]:
                    dms_data[column][site] = {}
                dms_data[column][site][mutant] = score
    return dms_data

def parse_gofasta_mutations(mutations_file):
    """Parse the gofasta CSV file to extract amino acid mutations
    
    Args:
        mutations_file (str): Path to CSV file containing mutations
        
    Returns:
        dict: Dictionary mapping sequence names to lists of mutation dictionaries
    """
    mutations_data = {}
    df = pd.read_csv(mutations_file)
    
    for _, row in df.iterrows():
        consensus = row['query']
        mutations_data[consensus] = []
        
        # Check if mutations column has a valid string value
        if pd.isna(row['mutations']) or not isinstance(row['mutations'], str):
            continue
            
        # Extract amino acid mutations
        aa_pattern = r'aa:[^:]+:([A-Z])(\d+)([A-Z])'
        all_mutations = row['mutations'].split('|')
        
        for mutation in all_mutations:
            if mutation.startswith('aa:'):
                match = re.search(aa_pattern, mutation)
                if match:
                    ref, pos, mutant = match.groups()
                    mutations_data[consensus].append({
                        'ref': ref,
                        'pos': int(pos),
                        'mutant': mutant
                    })
    
    return mutations_data

def calculate_dms_scores(dms_data, mutations_data):
    for sequence in mutations_data.keys():
        for idx, mutation in enumerate(mutations_data[sequence]):
            site = mutation['pos']
            mutant = mutation['mutant']
            mutations_data[sequence][idx]['dms_scores'] = {}
            # Get antibody_set from dms_data if available
            if 'antibody_set' in dms_data:
                mutations_data[sequence][idx]['antibody_set'] = dms_data['antibody_set'].get(site, {}).get(mutant, '')
            for column in dms_data.keys():
                if column != 'antibody_set':  # Skip antibody_set as it's handled separately
                    if site in dms_data[column]:
                        if mutant in dms_data[column][site]:
                            dms_score = dms_data[column][site][mutant]
                            mutations_data[sequence][idx]['dms_scores'][column] = dms_score
    return mutations_data

def clean_dms_scores(dms_scores):
    """Remove sequences with no mutations

    Args:
        dms_scores (dict): DMS scores for each sequence

    Returns:
        dict: DMS scores for each sequence with no mutations removed
    """
    return {k: v for k, v in dms_scores.items() if v}

def combine_dms_files(dms_file_1, dms_file_2):
    dms_df_1 = pd.read_csv(dms_file_1)
    dms_df_2 = pd.read_csv(dms_file_2)
    dms_df_2['antibody'] = dms_df_2['antibody'] + " sera escape"
    dms_df_2 = dms_df_2.pivot_table(index=["site", "wildtype", "mutant", "antibody_set"], 
                                   columns="antibody", 
                                   values="escape").reset_index()
    
    dms_df_1.set_index(['site', 'wildtype', 'mutant'], inplace=True)
    dms_df_2.set_index(['site', 'wildtype', 'mutant'], inplace=True)
    combined_df = pd.concat([dms_df_1, dms_df_2], axis=1)
    combined_df.reset_index(inplace=True)
    return combined_df

def write_dms_scores(dms_scores, output_file):
    print(f"Writing DMS scores to {output_file}")
    rows = []
    for consensus, mutations in dms_scores.items():
        for mutation in mutations:
            row = {
                "Consensus": consensus,
                "ref": mutation["ref"],
                "pos": mutation["pos"],
                "mutant": mutation["mutant"],
                "sra": consensus.split("_")[1] if len(consensus.split("_")) > 1 else "",
                "antibody_set": mutation.get("antibody_set", ""),
            }
            row.update(mutation["dms_scores"])
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dms-file", required=True, help="Input DMS data file")
    parser.add_argument("--dms-file-2", required=False, help="Input DMS data file")
    parser.add_argument('--h5_site_header', type=str, help='header of column with amino acid sites, in h5 sequential numbering')
    parser.add_argument('--mutation_file', type=str, help='CSV file with mutations generated by gofasta')
    parser.add_argument('--output_file', type=str, help='file to write DMS scores to')
    return parser.parse_args()

def main():
    args = parse_args()
    dms_df = combine_dms_files(args.dms_file, args.dms_file_2)

    # assert that each site has a unique mutant
    assert dms_df.groupby(['site', 'mutant']).size().max() == 1, "Each site must have a unique mutant"
    
    dms_data = read_dms_data(dms_df, h5_site_header=args.h5_site_header)
    mutations_data = parse_gofasta_mutations(args.mutation_file)
    dms_scores = calculate_dms_scores(dms_data, mutations_data)
    dms_scores = clean_dms_scores(dms_scores)
    write_dms_scores(dms_scores, args.output_file)

if __name__ == "__main__":
    main()