#########################################################################
# File Name: 1-get_pg_pair_FigR.r
# Author: Jun-hui Liu
# Last Modified: 2024-6-11
# Description: Run correlated Peak-Gene calling
# Usage: Fi-Cell
#########################################################################
# library(Seurat)
# library(FigR)
# library(SummarizedExperiment)
# library(foreach)
# library(doParallel)
library(getopt)
library(dplyr)

command<-matrix(c(
    'help', 'h', 0,'logical', 'help.',
    'cell', 'c', 1, 'character', 'cell type to run',
    'rna_seuret_obj', 'r', 1, 'character', 'seuret .rds file.',
    'atac_se_obj', 'a', 1, 'character', 'summarizedexperiment .rds file.',
    'windowPadSize','w',1,'numeric','base pairs padded on either side of gene TSS',
    'outdir','o',1,'character','outdir of correlated Peak-Gene pairs'
    ),byrow=T,ncol=5)
args<-getopt(command)
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

# cell = 'CD8+_T'
# rna_seuret_obj = '/mu03_ext2/liujunhui/Fi-Cell-main/example/rna_immune/CD8+_T.rds'
# atac_se_obj = '/mu03_ext2/liujunhui/Fi-Cell-main/example/atac_immune/CD8+_T.rds'
# windowPadSize=250000
# outdir='/mu03_ext2/liujunhui/Fi-Cell-main/FiCell-immune'

rna <- readRDS(rna_seuret_obj)
atac <- readRDS(atac_se_obj)
atac <- atac[,colData(atac)$DonorRace == 'White']
atac <- atac[grep("^chr([1-9]|1[0-9]|2[0-2]|X|Y)-", rownames(atac), value = TRUE), ]

############2. data check
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
all_value <- GetAssayData(rna, layer = 'counts')
if (all(all_value == floor(all_value))){
    rna <- NormalizeData(rna)
}
rm(all_value)

############3 Filter low quality cells and low expressed genes
### 3.1 Filter low quality cells
rna <- Seurat::AddMetaData(rna, metadata = GetAssayData(rna, assay = "RNA", layer = "counts") %>% colSums(), col.name = "nCount_RNA")
rna <- Seurat::AddMetaData(rna, metadata = GetAssayData(rna, assay = "RNA", layer = "counts") %>% apply(2, function(x) sum(x > 0)), col.name = "nFeature_RNA")
                                                                                                        
# min_genes_percentile <- 0.05
# min_counts_percentile <- 0.05
# min_atac_counts_percentile <- 0.05
                                                                                                        
# min_genes <- quantile(rna$nFeature_RNA, min_genes_percentile)
# min_counts <- quantile(rna$nCount_RNA, min_counts_percentile)
                                                                                                        
# cat(paste0("Filtering RNA cells: nFeature_RNA > ", round(min_genes), ", nCount_RNA > ", round(min_counts), "\n"))

min_genes <- 200
min_counts <- 500
rna_filt_cells <- colnames(rna)[rna$nFeature_RNA > min_genes & rna$nCount_RNA > min_counts]
cat(paste0("Filtering RNA cells: nFeature_RNA > ", round(min_genes), ", nCount_RNA > ", round(min_counts), "\n"))
# ATAC: filter cells based on total counts per cell
# atac_counts_per_cell <- Matrix::colSums(assays(atac)$counts)
# min_atac_counts <- quantile(atac_counts_per_cell, min_atac_counts_percentile)
# cat(paste0("Filtering ATAC cells: total counts > ", round(min_atac_counts), "\n"))
# atac_filt_cells <- colnames(atac)[atac_counts_per_cell > min_atac_counts]
# take intersection of good cells in RNA and ATAC
# common_cells <- intersect(rna_filt_cells, atac_filt_cells)
common_cells <- rna_filt_cells
                                                                                                        
cat(paste0("Keeping ", length(common_cells), " cells after filtering.\n"))
                                                                                                        
# subset RNA and ATAC to common cells
rna <- subset(rna, cells = common_cells)
atac <- atac[, common_cells]

### 3.2 Filter low expressed genes

rna_counts <- GetAssayData(rna, layer = "counts")
n_cells <- ncol(rna_counts)
gene_expressed_cells <- Matrix::rowSums(rna_counts > 0)

# keep genes expressed in >3% of cells
threshold_cells <- ceiling(0.03 * n_cells)
genes_to_keep <- names(gene_expressed_cells)[gene_expressed_cells > threshold_cells]

cat(paste0("Keeping ", length(genes_to_keep), " genes expressed in >3% of ", n_cells, " cells.\n"))

# subset RNA to selected genes
rna <- subset(rna, features = genes_to_keep)

# clean up memory
RNAmat <- GetAssayData(rna)                                                                                                     
rm(rna, rna_counts, gene_expressed_cells, genes_to_keep)
gc()

############4. Run correlated Peak-Gene calling

cisCor <- runGenePeakcorr(ATAC.se = atac,
                      RNAmat = RNAmat,
                      genome = "hg38",
                      keepMultiMappingPeaks = TRUE,
                      nCores = 16, 
                      normalizeATACmat=FALSE, #normalization has been done before
                      windowPadSize= windowPadSize,
                      p.cut= NULL)

cisCor$PG <- paste0(cisCor$PeakRanges, ":", cisCor$Gene)
write.table(cisCor, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.mmpeaks_kept.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

cisCor <- read.table(paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.mmpeaks_kept.tsv'), header = TRUE)
cisCor <- cisCor %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs==max(rObs))
cisCor$PG <- paste0(cisCor$PeakRanges, ":", cisCor$Gene)
write.table(cisCor, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

# cisCorfilt <- cisCor %>% filter(pvalZ <= 0.05)
# write.table(cisCorfilt, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCorfilt.pval5e-2.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

# cisCor$fdr <- p.adjust(cisCor$pvalZ, method = "fdr")
# cisCorfilt <- cisCor %>% filter(fdr <= 0.05)
# write.table(cisCorfilt, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCorfilt.fdr5e-2.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)

# cisCorfilt <- cisCor %>% filter(fdr <= 0.01)
# write.table(cisCorfilt, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCorfilt.fdr1e-2.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)


# for (cell in cells){
#     cisCor <- read.table(paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.mmpeaks_kept.tsv'), header = TRUE)
#     cisCor <- cisCor %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs==max(rObs))
#     cisCor$PG <- paste0(cisCor$PeakRanges, ":", cisCor$Gene)
#     write.table(cisCor, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCor.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
#     cisCorfilt <- cisCor %>% filter(pvalZ <= 0.05)
#     write.table(cisCorfilt, paste0(outdir,'/precomputation/cisCor/', cell,'.cisCorfilt.pval5e-2.tsv'), col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
# }
