# install.packages(c("BiocManager", "remotes", "qlcMatrix", "ggforce", "assertthat"))
# BiocManager::install(c("Signac", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", 
#                        "biovizBase", "RcisTarget", "org.Hs.eg.db",
#                        "TxDb.Hsapiens.UCSC.hg38.knownGene", "RCy3"))
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("satijalab/seurat-data")
# remotes::install_github("jiang-junyao/IReNA")
# remotes::install_github('satijalab/seurat-wrappers')

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SeuratDisk)
library(IReNA)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RCy3)
library(GenomicRanges)
library(org.Hs.eg.db)
library(annotate)
library(stringr)
library(dplyr)
library(visNetwork)
library(htmlwidgets)

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
                            refdata = reference$celltype.l2,
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

ggsave(plot = p1, filename = "MS4A1_1.jpg", dpi = "retina")

p2 <- CoveragePlot(object = pbmc,
                   region = "LYZ",
                   features = "LYZ",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 8000,
                   extend.downstream = 5000)

ggsave(plot = p2, filename = "LYZ_1.jpg", dpi = "retina")

patchwork::wrap_plots(p1, p2, ncol = 1)

# For diagonal integration "cca" analysis
reference <- NormalizeData(reference)
reference <- FindVariableFeatures(reference)
reference <- ScaleData(reference)
reference <- RunPCA(reference)
reference <- RunUMAP(reference, dims = 1:50)
DefaultAssay(pbmc) <- "SCT"

# transfer.anchors <- FindTransferAnchors(reference = reference, query = pbmc, reduction = "cca")
transfer.anchors <- readRDS("transfer.anchors.rds")

predictions <- TransferData(anchorset = transfer.anchors, 
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

ggsave(plot = p1, filename = "MS4A1_2.jpg", dpi = "retina")

p2 <- CoveragePlot(object = pbmc,
                   region = "LYZ",
                   features = "LYZ",
                   expression.assay = "SCT",
                   idents = idents.plot,
                   extend.upstream = 8000,
                   extend.downstream = 5000)

ggsave(plot = p2, filename = "LYZ_2.jpg", dpi = "retina")

patchwork::wrap_plots(p1, p2, ncol = 1)


# network analysis --------------------------------------------------------


# All peak should be linked to genes
# pbmc <- LinkPeaks(object = pbmc,
#                   peak.assay = "peaks",
#                   expression.assay = "SCT")

pbmc <- readRDS("pbmc_network.rds")

link_of_peaks = as.data.frame(Links(pbmc))

interval1 = str_split_fixed(link_of_peaks$peak, '-', 3)

interval1 = as.data.frame(interval1)

colnames(interval1) = c("Chromosome", "Start", "End")

intervals = GenomicRanges::makeGRangesFromDataFrame(interval1)

txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

genes = genes(txdb)

## Peak locations and overlapping genes
annotateIntervals <- function(intervals, txdb){
  
  stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))
  anno = genes(txdb)
  olaps = findOverlaps(intervals, anno)
  mcols(olaps)$gene_id = genes$gene_id[subjectHits(olaps)]
  intervals_factor = factor(queryHits(olaps), levels=seq_len(queryLength(olaps)))
  intervals$gene_id = splitAsList(mcols(olaps)$gene_id, intervals_factor)
  intervals
  
}  

myAnnotation <- as.data.frame(annotateIntervals(intervals, txdb))

myDf_master <- data.frame()
for (i in 1:length(myAnnotation$gene_id)) {
  
  # if the gene list is not empty...
  if(length(c(na.omit(myAnnotation$gene_id[i])[[1]])) != 0) {
    
    # annotate the interval and copy into a myDf data.frame
    myDf <- data.frame(myAnnotation$seqnames[i], 
                       myAnnotation$start[i], 
                       myAnnotation$end[i], 
                       toString(unname(getSYMBOL(c(na.omit(myAnnotation$gene_id[i])[[1]]), data='org.Hs.eg'))))
    
    # append tge myDF annotations with rbind into the myDf_master
    myDf_master <- rbind(myDf_master, myDf)
  }
}

myDf_header <- c("Chromosome", "Start", "End", "Genes")
names(myDf_master) <- myDf_header

myDf_master[c('Gene1','Gene_Other')] = str_split_fixed(myDf_master$Genes, ',', 2)

myDf_master$peaks = paste(myDf_master$Chromosome, 
                          myDf_master$Start, 
                          myDf_master$End, 
                          sep='-')


## Network Creation
link_of_genes = merge(x = link_of_peaks,
                      y = myDf_master, 
                      by.x = 'peak', 
                      by.y = "peaks")

link_only_trancription_factors = merge(x = link_of_genes,
                                       y = Tranfac201803_Hs_MotifTFsF,
                                       all.x = F, 
                                       all.y = F, 
                                       by.x = 'Gene1', 
                                       by.y = "TFs")


link_only_trancription_factors["interaction"]<-NA
for(i in 1:length(link_only_trancription_factors$score)){
  
  if(link_only_trancription_factors$score[i]<=0){
    link_only_trancription_factors$interaction[i]<-'inhibits'
  }else if (link_only_trancription_factors$score[i]>= 0){
    link_only_trancription_factors$interaction[i]<-'activates' 
  }else{
    link_only_trancription_factors$interaction[i]<-'interacts'
  }
  
}

nodes_dup = list(append(link_only_trancription_factors$Gene1, link_only_trancription_factors$gene))

un <- unlist(nodes_dup)
node_no_dup <- Map(`[`, nodes_dup, relist(!duplicated(un), skeleton = nodes_dup))

index = 1
groups = c()
for (i in 1:length((node_no_dup[[1]]))){
  if(node_no_dup[[1]][i] %in% Tranfac201803_Hs_MotifTFsF$TFs){
    groups[i]='TF'
  }
  else{
    groups[i]='Gene' 
  }
}

nodes = data.frame(id = node_no_dup,group=groups, 
                   label = node_no_dup)
colnames(nodes) = c('id', 'group', 'label')

edges = data.frame(from = link_only_trancription_factors$Gene1,
                   to =link_only_trancription_factors$gene, 
                   interaction = link_only_trancription_factors$interaction) 

edges = distinct(edges)

edge.color2 <- c()
for (i in edges$interaction) {
  if (i == "activates") {
    edge.color2 <- c(edge.color2, 'red')
  }
  else if (i == "inhibits") {
    edge.color2 <- c(edge.color2, 'blue')
  }
}

edges$color = edge.color2

graph = visNetwork(nodes, edges, height = "500px") %>%
  #visIgraphLayout(layout = "layout_in_circle") %>%
  visIgraphLayout() %>%
  visNodes(size = 10) %>%
  visGroups(group = "Gene", color = "#66AAAA", shape = "square", 
            shadow = list(enabled = TRUE)) %>%
  visGroups(group = "TF", color = "#FF9900", shape = "triangle") %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = T)

saveWidget(graph, file = "gene_regulatory_network.html")


# network analysis 2 ------------------------------------------------------


source("required_functions.R")

motif1 <- Tranfac201803_Hs_MotifTFsF

seurat_with_time = readRDS("seurat_with_time.rds")

load(system.file("extdata", "test_clustering.rda", package = "IReNA"))

expression_profile = test_clustering[,-1]

expression_profile <- na.omit(expression_profile)

# K-means clustering
clustering <- clustering_Kmeans(expression_profile, K1 = 4) #cluster number

plot_kmeans_pheatmap(clustering,ModuleColor1 = c('#67C7C1','#5BA6DA','#FFBF0F','#C067A9'))

# Add Ensembl ID as the first column of clustering results
Kmeans_clustering_ENS <- add_ENSID(clustering, Spec1 = 'Hs')

weightMat <- GENIE3(as.matrix(seurat_with_time@assays$RNA@data), nCores = 10)

weightMat <- getLinkList(weightMat)

regulation <- weightMat[weightMat[,3] > 0.0004, ]

regulatory_relationships <- add_regulation_type(Kmeans_clustering_ENS, regulation)

motif1 = Tranfac201803_Hs_MotifTFsF

motifTF <- c()
for (i in 1:nrow(motif1)) {
  TF <- strsplit(motif1[i,5],';')[[1]]
  motifTF <- c(motifTF,TF)
}

regulatory_relationships <- regulatory_relationships[regulatory_relationships[,1] %in% motifTF,]

gtf <- read.delim("/path/to/Homo_sapiens.GRCh38.108.chr.gtf", header=FALSE, comment.char="#")

gtf[,1] <- paste0('chr',gtf[,1])
gene_tss <- get_tss_region(gtf,rownames(Kmeans_clustering_ENS))

### Identify differentially expressed genes related motifs
motif1 <- motifs_select(Tranfac201803_Hs_MotifTFsF, rownames(Kmeans_clustering_ENS)) # Kmeans_clustering_ENS was obtained in part1

fimodir <- "/path/to/fimo"
outputdir1 <- '/path/to/out_dir/'
motifdir <- "/path/to/motif/"
refdir <- "/path/to/hg38.fa"

find_motifs_targetgenes2(gene_tss,
                         motif1, 
                         refdir,
                         fimodir,
                         outputdir1,
                         motifdir)

### run fimo_all script in shell
shell_code <- paste0('sh ', outputdir1, 'fimo/fimoall.sh')
system(shell_code, wait = TRUE)

motif2 <- Tranfac201803_Hs_MotifTFsF
outputdir <- paste0(outputdir1,'fimo/')
fimo_regulation <- generate_fimo_regulation(outputdir,motif2)

filtered_regulatory_relationships <- filter_regulation_fimo(fimo_regulation, 
                                                            regulatory_relationships)

TFs_list <- network_analysis(filtered_regulatory_relationships,
                             Kmeans_clustering_ENS,TFFDR1 = 10,
                             TFFDR2 = 10)

plot_tf_network(TFs_list)

enrichment_KEGG <- enrich_module(Kmeans_clustering_ENS, 
                                 org.Hs.eg.db,
                                 enrich.db = 'KEGG',
                                 organism = 'hsa',
                                 fun_num = 10, 
                                 use_internal_data = FALSE)

plot_intramodular_network(TFs_list,
                          enrichment_KEGG,
                          layout = 'circle')

