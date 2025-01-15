# Fi-Cell

- **F**unctional **i**nference of phenotype **C**ritical c**ell**;

## **Overview**

<p align="center">
    <img src="pic/figure-guthub.png" width="550"/>
</p>


## **Installation**

Fi-Cell can be downloaded by cloning this repository via the commands:

```shell
git clone https://github.com/LockenLiu/Fi-Cell.git
cd Fi-Cell 
```

## **Prerequisites**

### Anaconda or Miniconda distribution

In order to void any issues with software versioning, Fi-Cell utilises conda environments to automatically install all necessary dependencies. If conda is not present on your system, please install [anaconda](https://www.anaconda.com) or [miniconda](https://docs.conda.io/en/latest/miniconda.html) following the official instructions.

Once conda is ready, users can create an environment to install dependencies that Fi-Cell needs through the following commands:

```shell
conda install -c conda-forge mamba
mamba env create --file envs/FiCell.yaml
conda activate FiCell
```

Fi-Cell integrates [FigR](https://github.com/buenrostrolab/FigR) and [ldsc](https://github.com/bulik/ldsc), please install them before running Fi-Cell.

```R
devtools::install_github("buenrostrolab/FigR")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

```shell
cd your/dir/of/Fi-Cell
git clone https://github.com/bulik/ldsc.git
cd ldsc
conda activate
conda env create --file environment.yml
```

## **Getting Started** 

**Quick Run Fi-Cell**

```shell
conda activate FiCell
snakemake --use-conda --snakefile Fi-Cell.snakefile --configfile config.yml -j
```

We recommend running with `-j` as it will use all available cores. Specifying `-j 4` will use up to 4 cores. 
