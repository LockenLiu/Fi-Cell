outdir=$1
baselinemodeldir=$2
cell=$3
celllddir=$4
trait=$5
traitgwassumstatsdir=$6

mkdir -p ${outdir}/res

python ./ldsc/ldsc.py \
--h2 ${traitgwassumstatsdir} \
--ref-ld-chr ${baselinemodeldir}/baseline.,${celllddir}/${cell}. \
--frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--overlap-annot  \
--out ${outdir}/res/${cell}-${trait} \
--print-coefficients  \
--print-delete-vals  

# mkdir ${outdir}/res_condition

# conditionlddir=/mu03_ext2/liujunhui/Fi-Cell-main/FiCell-trait/precomputation/ldscore
# conditionprefix=Astrocytes
# python ./ldsc/ldsc.py \
# --h2 ${traitgwassumstatsdir} \
# --ref-ld-chr ${baselinemodeldir}/baseline.,${conditionlddir}/${conditionprefix}.,${celllddir}/${cell}. \
# --out ${outdir}/res_condition/${cell}-${trait} \
# --overlap-annot  \
# --frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
# --w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
# --print-coefficients  \
# --print-delete-vals  
