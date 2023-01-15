
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Integrated analysis for scRNA and scATAC seq data

## Installation

First, conda environment needs to be created.

    conda env create -f cell2cell_workshop_env.yml

Download required files.

    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
    wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi
    wget https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat

Download required tool.

    wget https://meme-suite.org/meme/meme-software/5.4.1/meme-5.4.1.tar.gz
    tar -xf meme-5.5.0.tar.gz
    cd meme-5.5.0
    ./configure
    make
    make install

## Usage

Activate the environment.

    conda activate cell2cell_workshop_env

Start a new R session by typing `R` and download required packages.

``` r
install.packages(c("BiocManager", "remotes", "qlcMatrix", "ggforce", "assertthat"))
BiocManager::install(c("Signac", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", 
                       "biovizBase", "RcisTarget", "org.Hs.eg.db",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "RCy3"))
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github("satijalab/seurat-data")
remotes::install_github("jiang-junyao/IReNA")
```

For peak calling, Macs2 path needs to be specified. It can be found by
typing `which macs2`.
