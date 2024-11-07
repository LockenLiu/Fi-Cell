cell=$1
dir=$2

dir=${dir}
mkdir ${dir}/precomputation/ldscore

for i in `seq 1 22`;
do
    python ./ldsc/make_annot.py \
    --bed-file ${dir}/precomputation/bed/${cell}.bed \
    --bimfile ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.bim \
    --annot-file ${dir}/precomputation/ldscore/${cell}.${i}.annot.gz
    #echo ${cell}:${i}/22
done