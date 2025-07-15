dir=$1
chr=$2
condition_prefix=$3

mkdir ${dir}/precomputation/condition_ldscore
bedtools multiinter -i ${dir}/precomputation/bed/*.bed | \
awk -v N=$(ls ${dir}/precomputation/bed/*.bed | wc -l) \
    -v halfN=$(( ( $(ls ${dir}/precomputation/bed/*.bed | wc -l) + 1 ) / 2 )) \
    '$4 >= halfN {print $1, $2, $3, 1}' OFS='\t' > \
    ${dir}/precomputation/condition_ldscore/${condition_prefix}.bed

python ./ldsc/make_annot.py \
--bed-file ${dir}/precomputation/condition_ldscore/${condition_prefix}.bed \
--bimfile ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
--annot-file ${dir}/precomputation/condition_ldscore/${condition_prefix}.${chr}.annot.gz


python ./ldsc/ldsc.py \
--l2 \
--bfile ./data/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
--ld-wind-cm 1 \
--annot ${dir}/precomputation/condition_ldscore/${condition_prefix}.${chr}.annot.gz \
--thin-annot \
--out ${dir}/precomputation/condition_ldscore/${condition_prefix}.${chr} \
--print-snps ./data/ldsc/print_snps.txt

echo combined:${chr}/22