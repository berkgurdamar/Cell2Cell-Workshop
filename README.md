
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

    wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
    gunzip Homo_sapiens.GRCh38.108.chr.gtf.gz

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz

    curl -o pbmc_network.rds -L 'https://drive.google.com/uc?export=download&confirm=yes&id=1jWZAA2l6ePa4JtFv2eViwCmSgINI72Lk'
    curl -o transfer.anchors.rds -L 'https://drive.google.com/uc?export=download&confirm=yes&id=1Y8TufsX4A6NynBksrLPJ0Mgvs1Y6dBrr'

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
remotes::install_github('satijalab/seurat-wrappers')
```

For peak calling, Macs2 path needs to be specified. It can be found by
typing `which macs2`.
