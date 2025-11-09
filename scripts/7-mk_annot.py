from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from pybedtools import BedTool
import gzip

parser = argparse.ArgumentParser(description='Create annotations from BED files for all chromosomes')
parser.add_argument('--cell', action="store",
                    dest="cell", type=str,
                    help='cell name')
parser.add_argument('--outdir', action="store",
                    dest="outdir", type=str,
                    help='The output path')
parser.add_argument('--chr', action="store",
                    dest="chr", type=str,
                    help='chr')

args = parser.parse_args()
bedname = args.cell
bedfile_path = args.outdir + "/bed"
annot_path = args.outdir + "/ldscore"
bimfile_path = './data/ldsc/1000G_EUR_Phase3_plink'
numchr = args.chr

def make_annot_files(bed_for_annot, bimfile, annot_file):
    print('making annot file')
    df_bim = pd.read_csv(bimfile,
            delim_whitespace=True, usecols = [0,1,2,3], names = ['CHR','SNP','CM','BP'])
    iter_bim = [['chr'+str(x1), x2, x2, 1] for (x1, x2) in np.array(df_bim[['CHR', 'BP']])]
    bimbed = BedTool(iter_bim)
    annotbed = bimbed.intersect(bed_for_annot, wb=True)
    bp = [x.start for x in annotbed]
    score = [float(x.fields[7]) for x in annotbed]
    df_int = pd.DataFrame({'BP': bp, 'ANNOT':score})
    df_annot = pd.merge(df_bim, df_int, how='left', on='BP')
    df_annot.fillna(0, inplace=True)
    temp = df_annot[['ANNOT']].astype(float)
    df_annot = pd.concat([df_bim.iloc[:,[0,3,1,2]], temp], axis = 1)
    thin_annot = temp
    if annot_file.endswith('.gz'):
        with gzip.open(annot_file, 'wb') as f:
            thin_annot.to_csv(f, sep = "\t", index = False)
    else:
        thin_annot.to_csv(annot_file, sep="\t", index=False)


# for numchr in range(1, 23, 1):
bimfile = bimfile_path + "/" + "1000G.EUR.QC." + str(numchr) + ".bim"
annot_file = annot_path + "/" + bedname + "." + str(numchr) + ".annot.gz"
bedfile = bedfile_path + "/" + bedname + ".bed"
bed_for_annot = BedTool(bedfile).sort()

make_annot_files(bed_for_annot, bimfile, annot_file)
print("We are at chrom : " + str(numchr))
