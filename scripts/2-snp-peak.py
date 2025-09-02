#########################################################################
# File Name: 2-snp-peak.py
# Author: LockenLiu
# Last Modified: 2025-8-20
#########################################################################

import pandas as pd
import os
from pybedtools import BedTool
import argparse

parser = argparse.ArgumentParser(description="find snp under peak")
parser.add_argument('-cd', '--cisCor_dir', help="cisCor directory", type=str, required=True)
parser.add_argument('-c', '--cell', help="cell to SNP", type=str, required=True)
parser.add_argument('-sc', '--snp_coordinate', help="snp_coordinate", type=str, required=True)

args = parser.parse_args()
cisCor_dir = args.cisCor_dir
cell = args.cell
snp_coordinate = args.snp_coordinate

# cisCor_dir = './FiCell-Signac/precomputation/cisCor/Astro'
# cell = 'Astro'
# snp_coordinate = 'data/hm3SNP-coordinate-hg38.tsv'

def find_snps_in_peak_bedtools(cisCor_bed, snp_bed):
    intersected = cisCor_bed.intersect(snp_bed, wa=True, wb=True)
    intersected_df = pd.read_csv(intersected.fn, sep='\t', header=None, 
                                 names=['chr', 'start', 'end', 'Peak', 'SNP_chr', 'SNP_start', 'SNP_end', 'SNP'])
    grouped_snps = intersected_df.groupby(['chr', 'start', 'end'])['SNP'].apply(lambda x: ', '.join(set(x))).reset_index()
    return grouped_snps

snp_cord = pd.read_csv(snp_coordinate, sep='\t')
snp_cord['START'] = pd.to_numeric(snp_cord['START'], errors='coerce').fillna(0).astype(int)
snp_cord['END'] = pd.to_numeric(snp_cord['END'], errors='coerce').fillna(0).astype(int)
snp_cord['CHR'] = snp_cord['CHR'].astype(str)  
snp_bed = BedTool.from_dataframe(snp_cord[['CHR', 'START', 'END', 'SNP']])

print(f'Processing cell type: {cell}')
try:
    cisCor = pd.read_csv(os.path.join(cisCor_dir, 'signac_peak_gene_links.tsv'), sep='\t')
    cisCor[['chr', 'start', 'end']] = cisCor['peak'].str.extract(r'chr(\d+|X|Y)-(\d+)-(\d+)')
    cisCor['start'] = pd.to_numeric(cisCor['start'], errors='coerce').fillna(0).astype(int)
    cisCor['end'] = pd.to_numeric(cisCor['end'], errors='coerce').fillna(0).astype(int)
    cisCor['chr'] = cisCor['chr'].astype(str)
    cisCor_bed = BedTool.from_dataframe(cisCor[['chr', 'start', 'end', 'peak']])
    snp_results = find_snps_in_peak_bedtools(cisCor_bed, snp_bed)
    snp_results['chr'] = snp_results['chr'].astype(str)
    cisCor = cisCor.merge(snp_results, on=['chr', 'start', 'end'], how='left')
    cisCor.rename(columns={'SNP': 'SNPs'}, inplace=True)
    output_file = os.path.join(cisCor_dir, f'{cell}.snpped.tsv')
    cisCor.to_csv(output_file, sep='\t', index=False)
    print(f'Successfully saved: {output_file}')
except FileNotFoundError:
    print(f"Error: File {'signac_peak_gene_links.tsv'} not found in {cisCor_dir}.")
except pd.errors.EmptyDataError:
    print(f"Error: File {'signac_peak_gene_links.tsv'} is empty or corrupted.")
except Exception as e:
    print(f"Error processing file {'signac_peak_gene_links.cisCor.tsv'}: {e}")
