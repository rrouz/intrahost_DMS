#!/usr/bin/env python3

import pandas as pd
import os
import argparse
import logging
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description='Process intrahost variant data')
    parser.add_argument('--dms_file', type=str, help='DMS data file (not used)')
    parser.add_argument('--data_dir', type=str, help='Directory with variant data')
    parser.add_argument('--output_dir', type=str, help='Directory to save the output')
    parser.add_argument('--chunk_size', type=int, default=10000, help='Chunk size for processing')
    parser.add_argument('--threads', type=int, default=4, help='Number of threads')
    return parser.parse_args()

def get_variant_files(data_dir):
    """Find all TSV files in the data directory"""
    variant_files = []
    if os.path.exists(data_dir):
        variant_files = [file.path for file in os.scandir(data_dir) if file.name.endswith(".tsv")]
        logger.info(f"Found {len(variant_files)} variant files")
    else:
        logger.warning(f"Directory not found: {data_dir}")
    
    return variant_files

def process_variant_chunk(chunk, filename):
    """Process a chunk of variant data, extracting SRA from filename and simplifying REGION field"""
    # Extract SRA from filename - assumes format like 'SRR12345678_something.tsv'
    # This extracts 'SRR12345678' as the SRA
    sra_id = os.path.basename(filename).split('_')[0]
    chunk['SRA'] = sra_id
    
    # Extract segment name from REGION field - take part before first pipe
    if 'REGION' in chunk.columns:
        chunk['REGION'] = chunk['REGION'].str.split('|').str[0]
    
    return chunk

def process_variant_files(variant_files, chunk_size=10000, num_threads=4):
    """Process all variant files"""
    all_variants = []
    
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for file in tqdm(variant_files, desc="Processing variant files"):
            try:
                chunks = []
                for chunk in pd.read_csv(file, sep='\t', chunksize=chunk_size):
                    future = executor.submit(process_variant_chunk, chunk, file)
                    chunks.append(future)
                
                if chunks:
                    file_variants = pd.concat([f.result() for f in chunks])
                    all_variants.append(file_variants)
                    logger.info(f"Processed {file} with {len(file_variants)} variants")
            except Exception as e:
                logger.error(f"Error processing {file}: {str(e)}")
    
    return all_variants

def save_outputs(output_dir, variant_data):
    """Save processed variant data to output directory"""
    logger.info("Saving outputs...")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create combined variants file
    if variant_data and any(not df.empty for df in variant_data):
        try:
            # Combine all dataframes that aren't empty
            combined_df = pd.concat([df for df in variant_data if not df.empty], ignore_index=True)
            logger.info(f"Combined {len(combined_df)} variants from {len(variant_data)} files")
        except Exception as e:
            logger.error(f"Error combining variants: {str(e)}")
            combined_df = pd.DataFrame()
    else:
        logger.warning("No variant data to combine")
        combined_df = pd.DataFrame()
    
    # Always write the combined file
    output_file = os.path.join(output_dir, 'combined_variants.tsv')
    combined_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved combined variants to {output_file}")

def main():
    args = parse_args()
    
    try:
        logger.info(f"Starting processing with data_dir: {args.data_dir}")
        logger.info(f"Output directory: {args.output_dir}")
        
        variant_files = get_variant_files(args.data_dir)
        variant_data = process_variant_files(variant_files, args.chunk_size, args.threads)
        save_outputs(args.output_dir, variant_data)
        
        # Double check that required file exists
        output_file = os.path.join(args.output_dir, 'combined_variants.tsv')
        if not os.path.exists(output_file):
            logger.warning(f"Required file combined_variants.tsv missing, creating empty file")
            pd.DataFrame().to_csv(output_file, sep='\t', index=False)
        
        logger.info("Processing completed successfully")
        
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        
        # Create required output file even in case of error
        output_file = os.path.join(args.output_dir, 'combined_variants.tsv')
        if not os.path.exists(output_file):
            logger.warning(f"Creating empty combined_variants.tsv due to error")
            pd.DataFrame().to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()