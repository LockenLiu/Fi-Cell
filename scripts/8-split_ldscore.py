import os
import sys
import re
import pandas as pd
#import pdb
import argparse


def split_M_files(file_ldscore, chromosome, prefix_genomic_annot, list_annotations):
	print("CHR={} | Writing .M and .M_5_50 files".format(chromosome))
	file_ldscore_base = re.sub(r"\.l2\.ldscore\.gz$", "", file_ldscore)
	file_M = "{}/{}.l2.M".format(cbmdir,file_ldscore_base)
	file_M_5_50 = "{}/{}.l2.M_5_50".format(cbmdir,file_ldscore_base)
	with open(file_M, "r") as fh_M, open(file_M_5_50, "r") as fh_M_5_50:
		list_M =  fh_M.readline().rstrip().split()
		list_M_5_50 = fh_M_5_50.readline().rstrip().split()
	assert(len(list_M) == len(list_M_5_50) == len(list_annotations))
	for i in range(len(list_annotations)):
		annotation = list_annotations[i]
		M = list_M[i]
		M_5_50 = list_M_5_50[i]
		file_out_M = "{}/{}.{}.l2.M".format(ctdir, annotation, chromosome)
		file_out_M_5_50 = "{}/{}.{}.l2.M_5_50".format(ctdir, annotation, chromosome)
		with open(file_out_M, "w") as fh_out_M, open(file_out_M_5_50, "w") as fh_out_M_5_50:
			fh_out_M.write(M)
			fh_out_M_5_50.write(M_5_50)
	print("CHR={} | DONE writing .M and .M_5_50 files".format(chromosome))

def split_ldscore_file_per_annotation(file_ldscore):
	print("Processing file_ldscore {}".format(file_ldscore))
	m = re.search(r"(.*)\.(\d{1,2})\.l2.ldscore.gz$", os.path.basename(file_ldscore))
	prefix_genomic_annot = m.groups()[0]
	chromosome = m.groups()[1]
	df = pd.read_csv("{}/{}".format(cbmdir,file_ldscore), sep="\t") # no index
	annotations_header = df.columns[3:].tolist() 
	annotations_clean = [re.sub(r"L2$", "", x) for x in annotations_header]
	split_M_files(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean)
	#split_annot_file(file_ldscore, chromosome, prefix_genomic_annot, annotations_clean)
	for counter, annotation in enumerate(annotations_header):
		annotation_clean = re.sub(r"L2$", "", annotation)
		file_out_ldscore = "{}/{}.{}.l2.ldscore.gz".format(ctdir, annotation_clean, chromosome)
		df[["CHR", "SNP", "BP", annotation]].to_csv(file_out_ldscore, sep="\t", index=False, compression="gzip")



###################################### MAIN ######################################

parser = argparse.ArgumentParser(description='Merging to get annot matrix')
parser.add_argument('-l', '--combined_ld_dir', help="combined_ld_dir", type=str, required=True)
parser.add_argument('-c', '--cell_ld_dir', help="cell_ld_dir", type=str, required=True)
parser.add_argument('-p', '--prefix', help="prefix", type=str, required=True)
parser.add_argument('-ch', '--chr', help="chr", type=str, required=True)

args = parser.parse_args()
cbmdir = args.combined_ld_dir
ctdir = args.cell_ld_dir
prefix = args.prefix
chro = args.chr

# cbmdir='FiCell-EXAMPLE/precomputation/Allcell-ldscore'
# ctdir='FiCell-EXAMPLE/precomputation/ldscore'
# prefix='Allcell'
# chro='1'

if not os.path.exists(ctdir):
    os.mkdir(ctdir)

one = prefix + '.' + str(chro) + '.l2.ldscore.gz'
split_ldscore_file_per_annotation(one)

checkfile = ctdir + '/' + prefix + '.' + str(chro) + '.splitedcheck.txt'
with open(checkfile, 'w') as f:
    print("done split {}".format(str(chro)))
    f.write('done split {}".format(str(chro))')