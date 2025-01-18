########################################################################################
################################### FUNCTIONS ##########################################
########################################################################################
include: "rules/functions.smk"

########################################################################################
################################### VARIABLES ##########################################
########################################################################################

# Where all the output will be saved
import os, re
BASE_OUTPUT_DIR = config['BASE_OUTPUT_DIR']
if not os.path.exists(BASE_OUTPUT_DIR):
    os.makedirs(BASE_OUTPUT_DIR)
CELL_TYPES = list(SC_SUMSTATS.keys())
CHROMOSOMES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
PREFIX = 'Allcell'
CONDITION_PREFIX = 'condition'

GWAS_TYPES = list(GWAS_SUMSTATS.keys())
# windowpadconsidered
windowpad_size = config['WINDOW_DEFINITION']['WINDOW_SIZE_KB']

BASELINE_MODEL_DIR = config['BASELINE_MODEL_DIR']

snp_coordinate_hg38 = config['SNP_COORDINATE']['hg38']
snp_coordinate_hg19 = config['SNP_COORDINATE']['hg19']

########################################################################################
################################### Target files ##########################################
########################################################################################

list_target_files = []
temp = expand("{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.tsv", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.snpped.tsv", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/SNPs/{cell_type}.cisCorfilt.pval5e-2.SNPs", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/bed/{cell_type}.bed", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.annot.gz", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, chromosome = CHROMOSOMES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.annot", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, PREFIX=PREFIX, chromosome = CHROMOSOMES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.l2.ldscore.gz", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, PREFIX=PREFIX, chromosome = CHROMOSOMES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.l2.ldscore.gz", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, chromosome = CHROMOSOMES)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/precomputation/condition_ldscore/{CONDITION_PREFIX}.{chromosome}.l2.ldscore.gz", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, chromosome=CHROMOSOMES, CONDITION_PREFIX=CONDITION_PREFIX)
list_target_files.append(temp)
temp = expand("{BASE_OUTPUT_DIR}/res/{cell_type}-{gwas_type}.results", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, gwas_type = GWAS_TYPES)
list_target_files.append(temp)

########################################################################################
################################### PIPELINE ##########################################
########################################################################################

rule all:
    input:
        list_target_files


rule get_pg_pair_FigR:
    input:
        lambda wildcards: SC_SUMSTATS[wildcards.cell_type]["rna_path"],
        lambda wildcards: SC_SUMSTATS[wildcards.cell_type]["atac_path"]
    output:
        "{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.tsv"
    params:
        cell = lambda wildcards: SC_SUMSTATS[wildcards.cell_type]['id'],
        rna_seuret_obj = lambda wildcards: SC_SUMSTATS[wildcards.cell_type]["rna_path"],
        atac_se_obj = lambda wildcards: SC_SUMSTATS[wildcards.cell_type]["atac_path"],
        windowPadSize = windowpad_size,
        outdir = config["BASE_OUTPUT_DIR"]
    shell:
        "Rscript scripts/1-get_pg_pair_FigR.r -c {params.cell} -r {params.rna_seuret_obj} -a {params.atac_se_obj} -w {params.windowPadSize} -o {params.outdir}"


rule snp_peak:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.tsv", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.snpped.tsv"
    params:
        cell = "{cell_type}", 
        cisCor_dir = "{BASE_OUTPUT_DIR}/precomputation/cisCor",
        snp_coord = snp_coordinate_hg38
    shell:
        "python scripts/2-snp-peak.py -cd {params.cisCor_dir} -c {params.cell} -sc {params.snp_coord}"


rule get_cts_peak_gene_pair_and_snps_under_peak:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/cisCor/{cell_type}.cisCor.snpped.tsv", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/SNPs/{cell_type}.cisCorfilt.pval5e-2.SNPs"
    params:
        cell = "{cell_type}", 
        cisCor_dir = "{BASE_OUTPUT_DIR}/precomputation/cisCor",
    shell:
        "Rscript scripts/3-get-cts-peak-gene-pair-and-snps-under-peak.r -c {params.cell} -d {params.cisCor_dir}"


rule SNP2bed:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/SNPs/{cell_type}.cisCorfilt.pval5e-2.SNPs", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/bed/{cell_type}.bed"
    params:
        cell = "{cell_type}", 
        snp_dir = "{BASE_OUTPUT_DIR}/precomputation/SNPs",
        snp_coord = snp_coordinate_hg19
    shell:
        "Rscript scripts/4-SNP2bed.r -c {params.cell} -d {params.snp_dir} -s {params.snp_coord}"

rule mk_annot:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/bed/{cell_type}.bed", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, chromosome = CHROMOSOMES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.annot.gz"
    conda:
        "ldsc"
    params:
        dir = "{BASE_OUTPUT_DIR}", 
        cell = "{cell_type}",
        chr = "{chromosome}", 
    shell:
        "sh scripts/5-mk_annot.sh {params.cell} {params.dir} {params.chr}"


rule merge_annot:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.annot.gz", BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, PREFIX=PREFIX, cell_type=CELL_TYPES, chromosome = CHROMOSOMES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.annot"
    params:
        dir = "{BASE_OUTPUT_DIR}/precomputation/ldscore", 
        prefix = "{PREFIX}"
    shell:
        "python scripts/6-combine_ant.py -p {params.prefix} -d {params.dir}"


rule mk_ldscore:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.annot",  BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, PREFIX=PREFIX, cell_type=CELL_TYPES, chromosome = CHROMOSOMES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.l2.ldscore.gz"
    conda:
        "ldsc"
    params:
        dir = "{BASE_OUTPUT_DIR}", 
        chr = "{chromosome}", 
        prefix = "{PREFIX}"
    shell:
        "sh scripts/7-mk_ldscore.sh {params.dir} {params.chr} {params.prefix}"


rule split_ldscore:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore/{PREFIX}.{chromosome}.l2.ldscore.gz", 
               BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, PREFIX=PREFIX, cell_type=CELL_TYPES, chromosome=CHROMOSOMES)
    output:
        "{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.l2.ldscore.gz"
    params:
        combined_ld_dir=lambda wildcards: f"{BASE_OUTPUT_DIR}/precomputation/{PREFIX}-ldscore",
        cell_ld_dir="{BASE_OUTPUT_DIR}/precomputation/ldscore",
        prefix=lambda wildcards: PREFIX, 
        chr = "{chromosome}"
    shell:
        "python scripts/8-split_ldscore.py -l {params.combined_ld_dir} -c {params.cell_ld_dir} -p {params.prefix} -ch {params.chr}"


rule mk_condition:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/bed/{cell_type}.bed", 
                BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, chromosome = CHROMOSOMES, CONDITION_PREFIX=CONDITION_PREFIX)
    output:
         "{BASE_OUTPUT_DIR}/precomputation/condition_ldscore/{CONDITION_PREFIX}.{chromosome}.annot.gz",
         "{BASE_OUTPUT_DIR}/precomputation/condition_ldscore/{CONDITION_PREFIX}.{chromosome}.l2.ldscore.gz"
    conda:
        "ldsc"
    params:
        outdir = lambda wildcards: BASE_OUTPUT_DIR,
        chr = "{chromosome}",
        cprefix = "{CONDITION_PREFIX}"
    shell:
        "sh scripts/9-mk_condition.sh {params.outdir} {params.chr} {params.cprefix}"


rule ldsc:
    input:
        expand("{BASE_OUTPUT_DIR}/precomputation/ldscore/{cell_type}.{chromosome}.l2.ldscore.gz", 
               BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, cell_type=CELL_TYPES, chromosome=CHROMOSOMES),
        expand("{BASE_OUTPUT_DIR}/precomputation/condition_ldscore/{CONDITION_PREFIX}.{chromosome}.l2.ldscore.gz", 
               BASE_OUTPUT_DIR=BASE_OUTPUT_DIR, CONDITION_PREFIX=CONDITION_PREFIX, chromosome=CHROMOSOMES),
        lambda wildcards: GWAS_SUMSTATS[wildcards.gwas_type]["path"]
    output:
         "{BASE_OUTPUT_DIR}/res/{cell_type}-{gwas_type}.results"
    conda:
        "ldsc"
    params:
        outdir=lambda wildcards: BASE_OUTPUT_DIR,
        baselinemodeldir=lambda wildcards: BASELINE_MODEL_DIR,
        cell=lambda wildcards: wildcards.cell_type,
        celllddir=lambda wildcards: f"{BASE_OUTPUT_DIR}/precomputation/ldscore",
        trait=lambda wildcards: wildcards.gwas_type,
        traitgwassumstatsdir=lambda wildcards: GWAS_SUMSTATS[wildcards.gwas_type]["path"],
        conditionlddir=lambda wildcards: f"{BASE_OUTPUT_DIR}/precomputation/condition_ldscore",
        conditionprefix=lambda wildcards: CONDITION_PREFIX
    shell:
        "sh scripts/10-ldsc-condition.sh {params.outdir} {params.baselinemodeldir} {params.cell} {params.celllddir} {params.trait} {params.traitgwassumstatsdir} {params.conditionlddir} {params.conditionprefix}"