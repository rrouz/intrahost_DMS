#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from tqdm import tqdm
from functools import lru_cache

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description='Process intrahost variant data')
    parser.add_argument('--dms_file', type=str, help='DMS data file')
    parser.add_argument('--data_dir', type=str, help='Directory with variant data')
    parser.add_argument('--output_dir', type=str, help='Directory to save the output')
    parser.add_argument('--chunk_size', type=int, default=10000, help='Chunk size for processing')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    return parser.parse_args()

def get_variant_files(data_dir):
    ha_variant_files, non_ha_variant_files = [], []
    for file in os.scandir(data_dir):
        if file.name.endswith(".tsv"):
            full_path = file.path
            if "HA" in file.name:
                ha_variant_files.append(full_path)
            else:
                non_ha_variant_files.append(full_path)
    
    logger.info(f"Found {len(ha_variant_files)} HA variant files and {len(non_ha_variant_files)} non-HA variant files")
    return ha_variant_files, non_ha_variant_files

def process_variant_chunk(chunk, filename):
    chunk['sra'] = os.path.basename(filename).split('_')[0]
    if 'REGION' in chunk.columns:
        chunk['REGION'] = chunk['REGION'].str.split('|').str[0]
    return chunk

def process_variant_files(variant_files, chunk_size=10000):
    variant_data = []
    with ThreadPoolExecutor() as executor:
        for file in variant_files:
            chunks = []
            for chunk in pd.read_csv(file, sep='\t', chunksize=chunk_size):
                future = executor.submit(process_variant_chunk, chunk, file)
                chunks.append(future)
            variant_data.append(pd.concat([f.result() for f in chunks]))
    return variant_data

@lru_cache(maxsize=1)
def read_dms_data(dms_file, h5_site_header='sequential_site'):
    logger.info(f"Reading DMS data from {dms_file}")
    dms_df = pd.read_csv(dms_file)
    
    columns_to_retrieve = [col for col in dms_df.columns 
                         if col not in ['Consensus', 'ref', 'mutant', 'sra', 'pos', 
                                       h5_site_header, 'nt changes to codon', 'antibody_set']]
    
    dms_data = {}
    for column in tqdm(columns_to_retrieve, desc="Processing DMS columns"):
        dms_data[column] = {}
        for _, row in dms_df.iterrows():
            if pd.isna(row[h5_site_header]) or pd.isna(row['mutant']):
                continue
                
            position = row[h5_site_header]
            mutation = row['mutant']
            score = row[column]
            
            if position not in dms_data[column]:
                dms_data[column][position] = {}
            
            dms_data[column][position][mutation] = score
    
    logger.info(f"Processed {len(columns_to_retrieve)} DMS score columns")
    return dms_data

def calculate_dms_scores_for_chunk(chunk, dms_data):
    result = chunk.copy()
    
    for column in dms_data.keys():
        scores = []
        for _, row in chunk.iterrows():
            try:
                scores.append(dms_data[column][row['POS_AA']][row['ALT_AA']])
            except (KeyError, TypeError):
                scores.append(np.nan)
        
        result[column] = scores
    
    return result

def calculate_dms_scores(dms_data, variant_data, num_threads=4, chunk_size=1000):
    logger.info("Calculating DMS scores...")
    processed_data = []
    
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for df in variant_data:
            chunks = np.array_split(df, max(1, len(df) // chunk_size))
            futures = [executor.submit(calculate_dms_scores_for_chunk, chunk, dms_data) 
                     for chunk in chunks]
            processed_chunks = [future.result() for future in futures]
            processed_data.append(pd.concat(processed_chunks))
    
    return processed_data

def save_outputs(output_dir, ha_variant_dms_data, non_ha_variant_data):
    logger.info("Saving outputs...")
    os.makedirs(output_dir, exist_ok=True)
    
    for idx, df in enumerate(ha_variant_dms_data):
        output_file = os.path.join(output_dir, f'ha_variant_{idx}.tsv')
        df.to_csv(output_file, sep='\t', index=False)
    
    for idx, df in enumerate(non_ha_variant_data):
        output_file = os.path.join(output_dir, f'non_ha_variant_{idx}.tsv')
        df.to_csv(output_file, sep='\t', index=False)
    
    combined_df = pd.concat(ha_variant_dms_data + non_ha_variant_data)
    combined_file = os.path.join(output_dir, 'combined_variants.tsv')
    combined_df.to_csv(combined_file, sep='\t', index=False)

def main():
    args = parse_args()
    
    try:
        ha_variant_files, non_ha_variant_files = get_variant_files(args.data_dir)
        
        ha_variant_data = process_variant_files(ha_variant_files, args.chunk_size)
        non_ha_variant_data = process_variant_files(non_ha_variant_files, args.chunk_size)
        
        dms_data = read_dms_data(args.dms_file)
        
        ha_variant_dms_data = calculate_dms_scores(
            dms_data, 
            ha_variant_data, 
            args.threads, 
            args.chunk_size
        )
        
        save_outputs(args.output_dir, ha_variant_dms_data, non_ha_variant_data)
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise

if __name__ == '__main__':
    main()