#########################################################################
# File Name: 3-get-cts-peak-gene-pair-and-snps-under-peak.r
# Author: Jun-hui Liu
# Last Modified: 2025-6-10
# Description: get cts peak gene pair and snps under peak
# Usage: Fi-Cell
#########################################################################
library(dplyr)
library(getopt)
library(tidyr)

command<-matrix(c(
    'help', 'h', 0,'logical', 'help.',
    'cell', 'c', 1, 'character', 'cell type to run',
    'cisCordir','d',1,'character','dir of cisCor'
    ),byrow=T,ncol=5)
args<-getopt(command)

if(!is.null(args$help)){
    help=getopt(command, usage = T)
    cat(paste(help,"get cts peak gene pair from cisCor pair.\n"))
}

cell=args$cell
dir=args$cisCordir
if (dir.exists(paste0(getwd(),"/",dir))){
  dir = paste0(getwd(),"/",dir)
}

# outdir='./FiCell-3_no_filt'
# dir='./FiCell-3_no_filt/precomputation/cisCor'
# cells = gsub('.rds', '', list.files('./example/rna_3/'))
# for (cell in cells){
#############1. 
cisCor <- read.table(paste0(dir, '/', cell, '.snpped.tsv'), header = TRUE, sep = '\t')
# sclinker <- read.table(paste0('/mu03_ext2/liujunhui/Fi-Cell-main/sclinker_immune/rna_GSE194122/', cell, '_Ltype.txt'))
# sclinker <- read.table(paste0('./sclinker/rna_GSE207334/', cell, '_Lcluster.txt'))
# names(sclinker) <- c('Gene', 'score')
# sclinker <- sclinker[sclinker$score > 0, ]
# cisCor <- na.omit(merge(cisCor, sclinker, by = 'Gene', all.x = TRUE))

#############2. 
cisCorfilt <- cisCor %>% filter(pvalue <= 0.05)
write.table(cisCorfilt, paste0(dir, '/', cell,'.pval5e-2.snpped.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
cisCorfilt_expanded <- cisCorfilt %>%
    separate_rows(SNPs, sep = ",\\s*") %>%
    filter(SNPs != "")
# snp_list <- cisCorfilt_expanded[,c('SNPs', 'score')]
snp_list <- data.frame(cisCorfilt_expanded$SNPs,1)
outdir <- paste0(dir, '/../SNPs/')
if (dir.exists(paste0(dir,"/../SNPs/")) == FALSE){
    outdir = paste0(dir,"/../SNPs/")
    dir.create(outdir)
}
write.table(snp_list, paste0(dir, '/../../SNPs/', cell, '.pval5e-2.SNPs'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#############3. 
cisCor$fdr <- p.adjust(cisCor$pvalue, method = "fdr")
cisCorfilt <- cisCor %>% filter(fdr <= 0.1)
write.table(cisCorfilt, paste0(dir, '/', cell,'.fdr1e-1.snpped.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
cisCorfilt_expanded <- cisCorfilt %>%
    separate_rows(SNPs, sep = ",\\s*") %>%
    filter(SNPs != "")
snp_list <- data.frame(cisCorfilt_expanded$SNPs,1)
# snp_list <- cisCorfilt_expanded[,c('SNPs', 'score')]
outdir <- paste0(dir, '/../SNPs/')
if (dir.exists(paste0(dir,"/../SNPs/")) == FALSE){
    outdir = paste0(dir,"/../SNPs/")
    dir.create(outdir)
}
write.table(snp_list, paste0(dir, '/../../SNPs/', cell, '.fdr1e-1.SNPs'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#############4. 
# cisCorfilt <- cisCor %>% filter(fdr <= 0.01)
# write.table(cisCorfilt, paste0(dir, '/', cell,'.cisCorfilt.fdr1e-2.snpped.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
# cisCorfilt_expanded <- cisCorfilt %>%
#     separate_rows(SNPs, sep = ",\\s*") %>%
#     filter(SNPs != "")
# snp_list <- data.frame(cisCorfilt_expanded$SNPs,1)
# outdir <- paste0(dir, '/../SNPs/')
# if (dir.exists(paste0(dir,"/../SNPs/")) == FALSE){
#     outdir = paste0(dir,"/../SNPs/")
#     dir.create(outdir)
# }
# write.table(snp_list, paste0(dir, '/../SNPs/', cell, '.cisCorfilt.fdr1e-2.SNPs'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

#################.1
# }