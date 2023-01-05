# install.packages(c("BiocManager", "hdf5r", "remotes", "qlcMatrix", "ggforce"))
# BiocManager::install(c("Signac", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "biovizBase"))
# remotes::install_github("mojaveazure/seurat-disk")

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)

set.seed(1234)

setwd("/path/to/files/")

# load the RNA and ATAC data
counts <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`,
                           assay = "RNA")

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,
                                       sep = c(":", "-"),
                                       fragments = fragpath,
                                       annotation = annotation)

DefaultAssay(pbmc) <- "ATAC"

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

VlnPlot(object = pbmc,
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 4,
        pt.size = 0)

# filter out low quality cells
pbmc <- subset(x = pbmc,
               subset = nCount_ATAC < 100000 &
                 nCount_RNA < 25000 &
                 nCount_ATAC > 1000 &
                 nCount_RNA > 1000 &
                 nucleosome_signal < 2 &
                 TSS.enrichment > 1)

# call peaks using MACS2
peaks <- CallPeaks(pbmc, outdir = "/desired/path", macs2.path = "/path/to/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(fragments = Fragments(pbmc),
                              features = peaks,
                              cells = colnames(pbmc))

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,
                                        fragments = fragpath,
                                        annotation = annotation)

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)

# load PBMC reference
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(reference = reference,
                                        query = pbmc,
                                        normalization.method = "SCT",
                                        reference.reduction = "spca",
                                        recompute.residuals = FALSE,
                                        dims = 1:50)

predictions <- TransferData(anchorset = transfer_anchors, 
                            efdata = reference$celltype.l2,
                            weight.reduction = pbmc[['pca']],
                            dims = 1:50)

pbmc <- AddMetaData(object = pbmc,
                    metadata = predictions)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                  "NK Proliferating", "gdT",
                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                  "CD14 Mono", "CD16 Mono",
                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(object = pbmc,
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:50, 2:40),
                                modality.weight.name = "RNA.weight",
                                verbose = TRUE)

# build a joint UMAP visualization
pbmc <- RunUMAP(object = pbmc,
                nn.name = "weighted.nn",
                assay = "RNA",
                verbose = TRUE)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes

pbmc <- LinkPeaks(object = pbmc,
                  peak.assay = "peaks",
                  expression.assay = "SCT",
                  genes.use = c("LYZ", "MS4A1"))

idents.plot <- c("B naive", "B intermediate", "B memory",
                 "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")

p1 <- CoveragePlot(object = pbmc,
                   region = "MS4A1",
                   features = "MS4A1",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 500,
                   extend.downstream = 10000)

p2 <- CoveragePlot(object = pbmc,
                   region = "LYZ",
                   features = "LYZ",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 8000,
                   extend.downstream = 5000)

patchwork::wrap_plots(p1, p2, ncol = 1)

# For diagonal integration "cca" analysis
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference)
reference <- RunUMAP(reference, dims = 1:50)
DefaultAssay(pbmc) <- "SCT"

transfer.anchors <- FindTransferAnchors(reference = reference, query = pbmc, reduction = "cca")

predictions <- TransferData(anchorset = transfer_anchors, 
                            refdata = reference$celltype.l2,
                            weight.reduction = pbmc[['pca']],
                            dims = 1:50)

pbmc <- AddMetaData(object = pbmc,
                    metadata = predictions)

Idents(pbmc) <- "predicted.id"

pbmc <- FindMultiModalNeighbors(object = pbmc,
                                reduction.list = list("pca", "lsi"), 
                                dims.list = list(1:50, 2:40),
                                modality.weight.name = "RNA.weight",
                                verbose = TRUE)

# build a joint UMAP visualization
pbmc <- RunUMAP(object = pbmc,
                nn.name = "weighted.nn",
                assay = "RNA",
                verbose = TRUE)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

DefaultAssay(pbmc) <- "peaks"

# first compute the GC content for each peak
pbmc <- RegionStats(pbmc, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
pbmc <- LinkPeaks(object = pbmc,
                  peak.assay = "peaks",
                  expression.assay = "SCT",
                  genes.use = c("LYZ", "MS4A1"))

idents.plot <- c("B naive", "B intermediate", "B memory",
                 "CD14 Mono", "CD16 Mono", "CD8 TEM", "CD8 Naive")

p1 <- CoveragePlot(object = pbmc,
                   region = "MS4A1",
                   features = "MS4A1",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 500,
                   extend.downstream = 10000)

p2 <- CoveragePlot(object = pbmc,
                   region = "LYZ",
                   features = "LYZ",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 8000,
                   extend.downstream = 5000)

patchwork::wrap_plots(p1, p2, ncol = 1)
