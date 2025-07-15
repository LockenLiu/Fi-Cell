cell=$1
dir=$2
chr=$3

dir=${dir}
mkdir ${dir}/precomputation/ldscore

python ./ldsc/make_annot.py \
--bed-file ${dir}/precomputation/bed/${cell}.bed \
--bimfile ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
--annot-file ${dir}/precomputation/ldscore/${cell}.${chr}.annot.gz