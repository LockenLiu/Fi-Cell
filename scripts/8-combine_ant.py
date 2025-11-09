#########################################################################
# File Name: 8-combine_ant.py
# Last Modified: 2025-11-09
#########################################################################


import os
import gzip
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='to construct annot matrix')
parser.add_argument('-d', '--annot_dir', help="annot filepath", type=str, required=True)
parser.add_argument('-p', '--prefix', help="prefix name", type=str, required=True)
args = parser.parse_args()
outdir = args.annot_dir
prefix = args.prefix


os.chdir(outdir)

if not os.path.exists('../{}-ldscore'.format(prefix)):
    os.mkdir('../{}-ldscore'.format(prefix))

files = [f for f in os.listdir('./') if f.endswith('.annot.gz')]

cell_types = sorted(set(f.split('.')[0] for f in files))
predefined_columns = sorted(cell_types)

chromosome_files = {}
for f in files:
    chr_number = f.split('.')[1]
    if chr_number not in chromosome_files:
        chromosome_files[chr_number] = []
    chromosome_files[chr_number].append(f)

for chr_number, files in chromosome_files.items():
    data_frames = []
    cell_names = []
    
    for f in files:
        cell_name = f.split('.')[0]
        cell_names.append(cell_name)
        
        with gzip.open(f, 'rt') as file:
            data = file.read().splitlines()

        data = [line for line in data if line != 'ANNOT']
        
        df = pd.DataFrame(data, columns=[cell_name])
        data_frames.append(df)
    
    combined_df = pd.concat(data_frames, axis=1)
    combined_df = combined_df[predefined_columns]
    output_filename = '../{}-ldscore/{}.{}.annot'.format(prefix, prefix, chr_number)
    combined_df.to_csv(output_filename, index=False, sep='\t')
    with open(output_filename, 'rb') as f_in:
        with gzip.open(f"{output_filename}.gz",'wb') as f_out:
            f_out.writelines(f_in)