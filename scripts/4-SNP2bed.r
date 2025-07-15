#########################################################################
# File Name: 4-SNP2bed.r
# Author: Jun-hui Liu
# Last Modified: 2024-6-11
# Description: SNP2bed.r
# Usage: Fi-Cell
#########################################################################
library(dplyr)
library(getopt)

command<-matrix(c(
    'help', 'h', 0,'logical', 'help.',
    'cell', 'c', 1, 'character', 'cell type to run',
    'dir','d',1,'character','dir of snp',
    'snp_coordinate', 's', 1, 'character', 'snp coordinate file'
    ),byrow=T,ncol=5)
args<-getopt(command)

if(!is.null(args$help)){
    help=getopt(command, usage = T)
    cat(paste(help,"get cts peak gene pair from cisCor pair.\n"))
}
cell=args$cell
snpdir=args$dir
if (dir.exists(paste0(getwd(),"/",snpdir))){
    dir = paste0(getwd(),"/",snpdir)
}

snp_coordinate=args$snp_coordinate

bim_bed <- read.table(snp_coordinate, header = TRUE,sep = '\t')
bim_bed$END <- bim_bed$END + 1
rownames(bim_bed) <- bim_bed$SNP
snp <- read.table(paste0(snpdir, '/', cell, '.pval5e-2.SNPs'))
snp <- cbind(bim_bed[snp$V1,], snp$V2)
snp <- snp %>%
    distinct(SNP, .keep_all = TRUE)
snp <- na.omit(snp)

snp <- snp[,c(1,3,4,5)]
names(snp)[4] <- 'score'
##不知道加不加score
# snp$score <- 1
snp$CHR <- paste0('chr',snp$CHR)
if (dir.exists(paste0(snpdir, "/../bed/")) == FALSE){
    outdir = paste0(snpdir, "/../SNPs/")
    dir.create(outdir)
}
write.table(snp, paste0(snpdir, '/../bed/',cell,".rawbed"), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

