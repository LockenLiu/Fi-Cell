#########################################################################
# File Name: 3-get-cts-peak-gene-pair-and-snps-under-peak.r
# Author: Jun-hui Liu
# Last Modified: 2024-6-11
# Description: get cts peak gene pair and snps under peak
# Usage: Fi-Cell
#########################################################################
library(dplyr)
library(getopt)

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

outdir=args$outdir
if (dir.exists(paste0(getwd(),"/",outdir))){
  outdir = paste0(getwd(),"/",outdir)
}

#############1. make all Peak-Gene dataframe
cell_list <- gsub('.cisCor.tsv', '', list.files(dir, pattern = ".cisCor.tsv$"))
all_pg <- c()
for (c in cell_list) {
    cisCor_cell <- read.table(paste0(dir, '/', c, '.cisCor.tsv'), header = TRUE)
    all_pg <- unique(c(all_pg, cisCor_cell$PG))
}
pg_df <- as.data.frame(matrix(0, nrow = length(all_pg), ncol = length(cell_list)))
rownames(pg_df) <- all_pg
colnames(pg_df) <- cell_list
for (c in cell_list) {
    cisCor_cell <- read.table(paste0(dir, '/', c, '.cisCor.tsv'), header = TRUE)
    pg_df[rownames(pg_df) %in% cisCor_cell$PG, c] <- 1
}
gene <- c()
for (i in 1:length(all_pg)){
    gene[i] <- strsplit(all_pg[i], ':')[[1]][3]
}
peak <- c()
for (i in 1:length(all_pg)){
    peak[i] <- paste0(strsplit(all_pg[i], ':')[[1]][1], ':', strsplit(all_pg[i], ':')[[1]][2])
}
pg_df$gene <- gene
pg_df$peak <- peak
pg_df$PG <- rownames(pg_df)


all_psnp <- read.table(paste0(dir, '/', cell,'.cisCor.snpped.tsv'), fill = TRUE, sep = '\t', header = TRUE)
all_psnp <- all_psnp[,c('PeakRanges','SNPs')]
for (c in cell_list) {
    psnp <- read.table(paste0(dir, '/', c,'.cisCor.snpped.tsv'), fill = TRUE, sep = '\t', header = TRUE)
    psnp <- psnp[,c('PeakRanges','SNPs')]
    all_psnp <- rbind(all_psnp, psnp)
    all_psnp <- all_psnp %>% distinct(PeakRanges, .keep_all = TRUE)   
}
pg_df <- merge(pg_df, all_psnp, by.x = 'peak', by.y = 'PeakRanges', all.x = TRUE)
rownames(pg_df) <- pg_df$PG

pg_df$sum <- rowSums(pg_df[,2:(length(cell_list) + 1)])
pg_df <- pg_df[pg_df$sum == 1,]


ctspair <- rownames(pg_df)[which(pg_df[cell] == 1)]
cisCor <- read.table(paste0(dir, '/', cell, '.cisCor.snpped.tsv'), header = TRUE, sep = '\t')
cisCorfilt <- cisCor[which(cisCor$pvalZ < 0.05),]
cisCorfilt <- cisCorfilt[cisCorfilt$PG %in% ctspair,]
write.table(cisCorfilt, paste0(dir, '/', cell,'.cisCorfilt.pval5e-2.cts.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
valid_snps <- cisCorfilt$SNPs[cisCorfilt$SNPs != ""]
snp_list <- unlist(strsplit(valid_snps, ", "))
snp_list <- unique(trimws(snp_list))
outdir <- paste0(dir, '/../SNPs/')
if (dir.exists(paste0(dir,"/../SNPs/")) == FALSE){
    outdir = paste0(dir,"/../SNPs/")
    dir.create(outdir)
}
write.table(snp_list, paste0(dir, '/../SNPs/', cell, '.cisCorfilt.pval5e-2.SNPs'), sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)