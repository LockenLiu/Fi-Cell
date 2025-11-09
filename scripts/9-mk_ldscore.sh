dir=$1
chr=$2
prefix=$3


python ./ldsc/ldsc.py \
--l2 \
--bfile ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
--ld-wind-cm 1 \
--annot ${dir}/precomputation/${prefix}-ldscore/${prefix}.${chr}.annot \
--thin-annot \
--out ${dir}/precomputation/${prefix}-ldscore/${prefix}.${chr} \
--print-snps ./data/ldsc/print_snps.txt
echo combined:${chr}/22
