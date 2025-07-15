outdir=$1
baselinemodeldir=$2
cell=$3
celllddir=$4
trait=$5
traitgwassumstatsdir=$6
conditionlddir=$7
conditionprefix=$8

mkdir -p ${outdir}/res_${conditionprefix}

# python ./ldsc/ldsc.py \
# --h2 /mu03_ext2/liujunhui/sumstats_formatted/cns_trait/Intelligence.sumstats.gz \
# --ref-ld-chr /mu03_ext2/liujunhui/Fi-Cell-main/data/ldsc/1000G_Phase3_baseline_v1.2_ldscores/baseline.,/mu03_ext2/liujunhui/Fi-Cell-main/FiCell-trait/precomputation/ldscore/Somatostatin-interneurons.,/mu03_ext2/liujunhui/Fi-Cell-main/FiCell-trait/precomputation/ldscore/LAMP5-interneurons. \
# --out ./LAMP5_AST_IQ \
# --overlap-annot  \
# --frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
# --w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# --print-coefficients  \
# --print-delete-vals  

# python ./ldsc/ldsc.py \
# --h2 ${traitgwassumstatsdir} \
# --ref-ld-chr ${baselinemodeldir}/baseline.,${conditionlddir}/${conditionprefix}.,${celllddir}/${cell}. \
# --out ${outdir}/res_${conditionprefix}/${cell}-${trait} \
# --overlap-annot  \
# --frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
# --w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# --print-coefficients  \
# --print-delete-vals  

python ./ldsc/ldsc.py \
--h2 ${traitgwassumstatsdir} \
--ref-ld-chr ${baselinemodeldir}/baseline.,${celllddir}/${cell}. \
--out ${outdir}/res_${conditionprefix}/${cell}-${trait} \
--overlap-annot  \
--frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  \
--print-delete-vals  
