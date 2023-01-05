
<!-- README.md is generated from README.Rmd. Please edit that file -->

Integrated analysis for scRNA and scATAC seq data

# Installation

First, conda environment needs to be created.

    conda env create -f cell2cell_workshop_env.yml

Downloading required files.

    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
    wget https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat

# Usage

Activate the environment.

    conda activate cell2cell_workshop_env

Start a new R session by typing `R` and download required packages.

``` r
install.packages(c("BiocManager", "remotes", "qlcMatrix", "ggforce"))
BiocManager::install(c("Signac", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "biovizBase"))
remotes::install_github("mojaveazure/seurat-disk")
```

For peak calling, Macs2 path needs to be specified. It can be found by
typing `which macs2`.