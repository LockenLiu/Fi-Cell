#########################################################################
# File Name: 1-1-get_pg_pair_Signac.r
# Author: LockenLiu
# Last Modified: 2025-8-20
# This script is modified based on:  https://github.com/elizabethdorans/E2G_Method_Tutorials/blob/main/Signac/run_Signac_single_chromosome.R
#########################################################################

suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(data.table))

set.seed(1234)

parser <- ArgumentParser()

parser$add_argument("chromosome",
    help="[REQUIRED] chromosome whose genes are being tested, i.e. 1")
parser$add_argument("--seurat_object",
    help="[REQUIRED] Path to Seurat object preprocessed for Signac peak-gene linking (see Signac_preprocessing.R) [.rds]")
parser$add_argument("--signac_output_dir", default = ".",
    help = "Path to directory for output files")
parser$add_argument("--max_peak_TSS_distance", default = 5e+05,
    help = "Maximum distance (in base pairs) between TSS and peak")

args <- parser$parse_args()

chromosome = args$chromosome
seurat_object = args$seurat_object
signac_output_dir = args$signac_output_dir
max_peak_TSS_distance = as.numeric(args$max_peak_TSS_distance)

# Create output directory if needed
print(sprintf("Output directory: %s", signac_output_dir))
if (!dir.exists(signac_output_dir)) {
    dir.create(signac_output_dir, recursive = TRUE)
}

# Create outfile name
links_outfile <- sprintf("%s/chr%s.tsv", signac_output_dir, chromosome)

# Read in Seurat object
print("Reading in Seurat object!")
data <- readRDS(seurat_object)

# Define set of genes on chromosome for peak-gene linking
annot <- Annotation(object = data[["ATAC"]])
gene_coords = as.data.table(annot)[, c("seqnames", "gene_name")]
gene_coords = unique(gene_coords)
gene_set = gene_coords[seqnames == sprintf("chr%s", chromosome), ]$gene_name

sprintf("Linking across %s possible genes on chromosome %s!", length(gene_set), chromosome)

DefaultAssay(data) <- "ATAC"

data <- LinkPeaks(
    object = data,
    peak.assay = "ATAC",
    expression.assay = "SCT",
    genes.use = gene_set,
    distance = max_peak_TSS_distance,
    pvalue_cutoff = 0.05,
    score_cutoff = 0
)

out_df = as.data.frame(Links(data))[,c("peak", "gene", "score", "pvalue")]
write.table(out_df, links_outfile, sep = "\t", quote = FALSE, row.names = FALSE)
sprintf("Output to %s!", links_outfile)