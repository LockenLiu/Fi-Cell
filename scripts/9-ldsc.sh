outdir=$1
baselinemodeldir=$2
cell=$3
celllddir=$4
trait=$5
traitgwassumstatsdir=$6

mkdir ${outdir}/res

python ~/Biosofts/ldsc/ldsc.py \
--h2 ${traitgwassumstatsdir} \
--ref-ld-chr ${baselinemodeldir}/baseline.,${celllddir}/${cell}. \
--out ${outdir}/res/${cell}-${trait} \
--overlap-annot  \
--frqfile-chr ./data/ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
--w-ld-chr ./data/ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients  \
--print-delete-vals  