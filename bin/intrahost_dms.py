#!/usr/bin/env python3

import pandas as pd
import os
from tqdm import tqdm
import argparse
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict
from functools import lru_cache
import logging

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

def get_variant_files(data_dir: str) -> tuple[List[str], List[str]]:
    ha_variant_files, non_ha_variant_files = [], []
    for file in os.scandir(data_dir):  # Using scandir instead of listdir
        if file.name.endswith(".tsv"):
            full_path = file.path
            if "HA" in file.name:
                ha_variant_files.append(full_path)
            else:
                non_ha_variant_files.append(full_path)
    return ha_variant_files, non_ha_variant_files

def process_variant_chunk(chunk: pd.DataFrame, filename: str) -> pd.DataFrame:
    chunk['sra'] = filename.split('/')[-1].split('_')[0]
    return chunk

def process_variant_files(data_dir: str, variant_files: List[str], chunk_size: int = 10000) -> List[pd.DataFrame]:
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
def read_dms_data(dms_file: str, h5_site_header: str) -> Dict:
    logger.info("Reading DMS data...")
    dms_df = pd.read_csv(dms_file)
    columns_to_retrieve = dms_df.columns[3:].tolist()
    
    dms_data = {}
    for column in tqdm(columns_to_retrieve, desc="Processing DMS columns"):
        grouped = dms_df.groupby([h5_site_header, 'mutant'])[column].first()
        dms_data[column] = grouped.unstack().to_dict()
    
    return dms_data

def calculate_dms_scores_for_chunk(args: tuple) -> pd.DataFrame:
    chunk, dms_data = args
    result = chunk.copy()
    for column in dms_data.keys():
        scores = []
        for _, row in chunk.iterrows():
            try:
                scores.append(dms_data[column][row['POS_AA']][row['ALT_AA']])
            except KeyError:
                scores.append(np.nan)
        result[column] = scores
    return result

def calculate_dms_scores(dms_data: Dict, variant_data: List[pd.DataFrame], 
                        num_threads: int = 4, chunk_size: int = 1000) -> List[pd.DataFrame]:
    logger.info("Calculating DMS scores...")
    processed_data = []
    
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for df in variant_data:
            chunks = np.array_split(df, max(1, len(df) // chunk_size))
            futures = [executor.submit(calculate_dms_scores_for_chunk, (chunk, dms_data)) 
                      for chunk in chunks]
            processed_chunks = [future.result() for future in futures]
            processed_data.append(pd.concat(processed_chunks))
    
    return processed_data

def save_outputs(output_dir: str, ha_variant_dms_data: List[pd.DataFrame], 
                non_ha_variant_data: List[pd.DataFrame]):
    logger.info("Saving outputs...")
    os.makedirs(output_dir, exist_ok=True)
    
    # Use parallel processing for file saving
    with ThreadPoolExecutor() as executor:
        # Save HA variant files
        futures = []
        for idx, variant_df in enumerate(ha_variant_dms_data):
            future = executor.submit(variant_df.to_csv, 
                                  os.path.join(output_dir, f'ha_variant_{idx}.tsv'),
                                  sep='\t', index=False)
            futures.append(future)
        
        # Save non-HA variant files
        for idx, variant_df in enumerate(non_ha_variant_data):
            future = executor.submit(variant_df.to_csv,
                                  os.path.join(output_dir, f'non_ha_variant_{idx}.tsv'),
                                  sep='\t', index=False)
            futures.append(future)
        
        # Wait for all files to be saved
        for future in futures:
            future.result()
    
    # Save combined data
    combined_df = pd.concat(ha_variant_dms_data + non_ha_variant_data)
    combined_df.to_csv(os.path.join(output_dir, 'combined_variants.tsv'), 
                      sep='\t', index=False)

def main():
    args = parse_args()
    
    try:
        ha_variant_files, non_ha_variant_files = get_variant_files(args.data_dir)
        
        ha_variant_data = process_variant_files(args.data_dir, ha_variant_files, 
                                              args.chunk_size)
        non_ha_variant_data = process_variant_files(args.data_dir, non_ha_variant_files, 
                                                  args.chunk_size)
        
        dms_data = read_dms_data(args.dms_file, 'sequential_site')
        
        ha_variant_dms_data = calculate_dms_scores(dms_data, ha_variant_data, 
                                                 args.threads, args.chunk_size)
        
        save_outputs(args.output_dir, ha_variant_dms_data, non_ha_variant_data)
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}")
        raise

if __name__ == '__main__':
    main()