#########################################################################
# File Name: 1-get_pg_pair_FigR.r
# Author: Jun-hui Liu
# Last Modified: 2024-6-11
# Description: Run correlated Peak-Gene calling
# Usage: Fi-Cell
#########################################################################
library(Seurat)
library(FigR)
library(SummarizedExperiment)
library(foreach)
library(doParallel)
library(getopt)

command<-matrix(c(
    'help', 'h', 0,'logical', 'help.',
    'cell', 'c', 1, 'character', 'cell type to run',
    'rna_seuret_obj', 'r', 1, 'character', 'seuret .rds file.',
    'atac_se_obj', 'a', 1, 'character', 'summarizedexperiment .rds file.',
    'windowPadSize','w',1,'numeric','base pairs padded on either side of gene TSS',
    'outdir','o',1,'character','outdir of correlated Peak-Gene pairs'
    ),byrow=T,ncol=5)
args<-getopt(command)
clump_file = args$clump_file
if(!is.null(args$help) ){
    help=getopt(command, usage = T)
    cat(paste(help,"calling peak-gene pair in cell types.\n"))
}

#############1. load data
cell = args$cell
rna_seuret_obj = args$rna_seuret_obj
if (file.exists(paste0(getwd(),"/",rna_seuret_obj))){
  rna_seuret_obj = paste0(getwd(),"/",rna_seuret_obj)
}
atac_se_obj = args$atac_se_obj
if (file.exists(paste0(getwd(),"/",atac_se_obj))){
  atac_se_obj = paste0(getwd(),"/",atac_se_obj)
}
windowPadSize=args$windowPadSize
outdir=args$outdir
if (dir.exists(paste0(getwd(),"/",outdir))){
  outdir = paste0(getwd(),"/",outdir)
}
rna <- readRDS(rna_seuret_obj)
atac <- readRDS(atac_se_obj)

#############2. data check
if (dim(rna)[[2]] != dim(atac)[[2]]){
    stop("cell counts in RNA and ATAC datasets are not equal, please embed two datasets and unify the colnames first.")
}
if (FALSE %in% (dim(rna)[[2]] == dim(atac)[[2]])){
    stop("colnames in RNA and ATAC datasets are consistent , please embed two datasets and unify the colnames first.")
}
# perform ATACmat and RNAmat normalization
all_value <- assays(atac)$counts@x
if (all(all_value == floor(all_value))){
    atac_seurat <- Seurat::CreateSeuratObject(counts = assays(atac)$counts)
    atac_seurat <- NormalizeData(atac_seurat)
    assays(atac)$counts <- GetAssayData(atac_seurat)
}
all_value <- GetAssayData(rna)
if (all(all_value == floor(all_value))){
    rna <- NormalizeData(rna)
}
rm(all_value)

#############3. Run correlated Peak-Gene calling
RNAmat <- GetAssayData(rna)
cisCor <- runGenePeakcorr(ATAC.se = atac,
                      RNAmat = RNAmat,
                      genome = "hg38",
                      keepMultiMappingPeaks = TRUE, #remove MultiMapped Peaks later
                      nCores = 32, 
                      normalizeATACmat=FALSE, #normalization has been done before
                      windowPadSize= 250000,
                      p.cut= NULL)

cisCor <- cisCor %>% dplyr::group_by(Peak) %>% dplyr::filter(pvalZ==min(pvalZ))
cisCor$PG <- paste0(cisCor$PeakRanges, ":", cisCor$Gene)
write.table(cisCor, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
cisCorfilt <- cisCor %>% filter(pvalZ <= 0.05)
write.table(cisCorfilt, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCorfilt.pval5e-2.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
