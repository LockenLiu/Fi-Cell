#########################################################################
# File Name: 4-get-sig-peak-gene-pair-and-snps-under-peak.r
# Last Modified: 2025-11-09
#########################################################################

library(dplyr)
library(getopt)
library(tidyr)

command<-matrix(c(
    'help', 'h', 0,'logical', 'help.',
    'cell', 'c', 1, 'character', 'cell type to run',
    'cts', 's', 1, 'character', 'gene_cts_score',
    'cisCordir','d',1,'character','dir of cisCor'
    ),byrow=T,ncol=5)
args<-getopt(command)

if(!is.null(args$help)){
    help=getopt(command, usage = T)
    cat(paste(help,"get cts peak gene pair from cisCor pairs.\n"))
}

cell=args$cell
cts_file=args$cts
dir=args$cisCordir
if (dir.exists(paste0(getwd(),"/",dir))){
  dir = paste0(getwd(),"/",dir)
}


#############1. read in cisCor and  cell-type-specificity file
cisCor <- read.table(paste0(dir, '/', cell, '.snpped.tsv'), header = TRUE, sep = '\t')
cts_score <- read.table(cts_file)
names(cts_score) <- c('Gene', 'score')
cts_score <- cts_score[cts_score$score > 0, ]
cisCor <- na.omit(merge(cisCor, cts_score, by.x = 'gene', by.y = 'Gene', all.x = TRUE))

#############2. 
cisCorfilt <- cisCor %>% filter(pvalue <= 0.05)
write.table(cisCorfilt, paste0(dir, '/', cell,'.pval5e-2.snpped.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
cisCorfilt_expanded <- cisCorfilt %>%
    separate_rows(SNPs, sep = ",\\s*") %>%
    filter(SNPs != "")
snp_list <- cisCorfilt_expanded[,c('SNPs', 'score')]
# snp_list <- data.frame(cisCorfilt_expanded$SNPs,1)
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
# snp_list <- data.frame(cisCorfilt_expanded$SNPs,1)
snp_list <- cisCorfilt_expanded[,c('SNPs', 'score')]
outdir <- paste0(dir, '/../SNPs/')
if (dir.exists(paste0(dir,"/../SNPs/")) == FALSE){
    outdir = paste0(dir,"/../SNPs/")
    dir.create(outdir)
}
write.table(snp_list, paste0(dir, '/../../SNPs/', cell, '.fdr1e-1.SNPs'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
