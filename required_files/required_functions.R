getfasta <- function(merged_footprints, fastadir) {
  fasta <- Biostrings::readBStringSet(fastadir, format = "fasta", nrec = -1L,
                                      skip = 0L, seek.first.rec = FALSE,
                                      use.names = TRUE)
  fasta1 <- c()
  for (i in 1:nrow(merged_footprints)) {
    name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                   "-", merged_footprints[i, 3])
    if (merged_footprints[i, 3] > length(fasta[[merged_footprints[i, 1]]])) {
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):length(fasta[[merged_footprints[i, 1]]])]))
      name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                     "-", length(fasta[[merged_footprints[i, 1]]]))
    } else{
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):merged_footprints[i, 3]]))
    }
    fasta1 <- c(fasta1, name, sequence)
  }
  fasta1 <- as.data.frame(fasta1)
  return(fasta1)
}

getfasta2 <- function(merged_footprints, fasta) {
  fasta1 <- c()
  for (i in 1:nrow(merged_footprints)) {
    name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                   "-", merged_footprints[i, 3])
    if (merged_footprints[i, 3] > length(fasta[[merged_footprints[i, 1]]])) {
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):length(fasta[[merged_footprints[i, 1]]])]))
      name <- paste0(">", merged_footprints[i, 1], ":", merged_footprints[i, 2],
                     "-", length(fasta[[merged_footprints[i, 1]]]))
    } else{
      sequence <- toupper(as.character(fasta[[merged_footprints[i, 1]]][
        (merged_footprints[i, 2] + 1):merged_footprints[i, 3]]))
    }
    fasta1 <- c(fasta1, name, sequence)
  }
  fasta1 <- as.data.frame(fasta1)
  return(fasta1)
}







get_merged_fasta2 = function (footprints, fastadir, distance = 4) 
{
  #validInput(footprints, "footprints", "df")
  #validInput(fastadir, "fastadir", "direxists")
  #validInput(distance, "distance", "numeric")
  merged_footprints <- merge_footprints(footprints, ditance = distance)
  fasta <- getfasta(merged_footprints, fastadir = fastadir)
  return(fasta)
}


find_motifs_targetgenes2 = function (gene_tss, motif, refdir, fimodir, outputdir1, Motifdir, 
                                     sequencedir = NULL, select_motif = T, use_nohup = F) 
{
  #validInput(refdir, "refdir", "fileexists")
  #validInput(Motifdir, "Motifdir", "direxists")
  #validInput(outputdir1, "outputdir", "direxists")
  if (str_ends(outputdir1, "/") == FALSE) {
    warning("the last character of outputdir1 is not \"/\"")
  }
  fasta <- Biostrings::readBStringSet(refdir, format = "fasta", 
                                      nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
  gene_tss[, 3] <- as.integer(gene_tss[, 3])
  gene_tss[, 4] <- as.integer(gene_tss[, 4])
  if (select_motif == T) {
    motif1 <- motifs_select(motif, gene_tss[, 1])
  }
  else {
    motif1 = motif
  }
  outputdir <- paste0(outputdir1, "fimo/")
  if (is.null(sequencedir)) {
    sequencedir <- paste0(outputdir1, "fasta/")
  }
  fimoall <- c()
  if (!"fasta" %in% dir(outputdir1)) {
    dir.create(paste0(outputdir1, "fasta"))
  }
  if (!"fimo" %in% dir(outputdir1)) {
    dir.create(paste0(outputdir1, "fimo"))
  }
  for (i in 1:nrow(gene_tss)) {
    if (use_nohup == TRUE) {
      fimo1 <- paste0("nohup sh ", outputdir1, "fimo/", 
                      gene_tss[i, 1], "/", "Fimo1.sh &")
    }
    else if (use_nohup == FALSE) {
      fimo1 <- paste0("sh ", outputdir1, "fimo/", gene_tss[i, 
                                                           1], "/", "Fimo1.sh ;")
    }
    else {
      stop("parameter use_nohup should be TRUE or FALSE")
    }
    fimoall <- c(fimoall, fimo1)
    dir.create(paste0(outputdir1, "fimo/", gene_tss[i, 1]))
    fasta1 <- getfasta2(gene_tss[i, 2:4], fasta)
    write.table(fasta1, paste0(outputdir1, "fasta/", gene_tss[i, 
                                                              1], ".fa"), col.names = F, row.names = F, quote = F)
    fimodir1 <- fimodir
    outputdir12 <- paste0(outputdir1, "fimo/", gene_tss[i, 
                                                        1], "/")
    outputdir11 <- paste0(outputdir, gene_tss[i, 1], "/")
    sequencedir1 <- paste0(sequencedir, gene_tss[i, 1], ".fa")
    find_motifs(motif1, step = 500, fimodir1, outputdir12, 
                outputdir11, Motifdir, sequencedir1)
  }
  fimoall <- as.data.frame(fimoall)
  write.table(fimoall, paste0(outputdir1, "fimo/", "fimoall.sh"), 
              quote = F, row.names = F, col.names = F)
}




plot_tf_network2 = function (TFs_list, layout = "grid", group.cols = NULL, title.name = NULL, 
                             vertex.size = 13, vertex.size.add = 3, vertex.label.color = "black", 
                             edge.label.color = "black", legend = TRUE, vertex.label.cex = 0.8, 
                             vertex.label.family = "ArialMT", frame.color = "white", arrow.size = 0.2, 
                             arrow.width = 0.5, edge.width = 1.8, edge.curved = 0, edge.color = c("#FDD1B0", 
                                                                                                  "#B3B3B3")) 
{
  #validInput(TFs_list, "TFs_list", "list")
  #validInput(vertex.size, "vertex.size", "numeric")
  #validInput(vertex.size.add, "vertex.size.add", "numeric")
  #validInput(legend, "legend", "logical")
  #validInput(vertex.label.cex, "vertex.label.cex", "numeric")
  #validInput(arrow.size, "arrow.size", "numeric")
  #validInput(arrow.width, "arrow.width", "numeric")
  #validInput(edge.width, "edge.width", "numeric")
  #validInput(edge.curved, "edge.curved", "numeric")
  if (length(edge.color) != 2) {
    stop("You should input two colors for edge.color")
  }
  if (is.null(group.cols)) {
    col1 <- c("#67C7C1", "#5BA6DA", "#FFBF0F", "#C067A9", 
              "#EF9951", "#E56145", "#C0C130", "#67C1E3", "#EF9951", 
              "#00BFC4", "#AEC7E8", "#E56145", "#2F4F4F")
  }
  else {
    col1 <- group.cols
  }
  network <- TFs_list[["FOSF_RegMTF_Cor_EnTFs"]]
  network$TFGroup <- as.integer(network$TFGroup)
  network$TargetGroup <- as.integer(network$TargetGroup)
  tfs <- network[, c("TFSymbol", "TFGroup")]
  target <- network[, c("TargetSymbol", "TargetGroup")]
  colnames(target) <- c("TFSymbol", "TFGroup")
  nodes <- rbind(tfs, target)
  edges <- network[, c("TFSymbol", "TargetSymbol", "Regulation", 
                       "Correlation")]
  edge_type <- network$Regulation
  colnames(nodes) <- c("name", "type")
  nodes <- nodes[!duplicated(nodes$name), ]
  colnames(edges) <- c("from", "to", "type", "weight")
  g <- igraph::graph_from_data_frame(edges, vertices = nodes, 
                                     directed = TRUE)
  tf_degree <- as.data.frame(table(network$TFSymbol))
  cut1 <- as.numeric(summary(as.data.frame(table(network$TFSymbol))[, 
                                                                    2])[2])
  cut2 <- as.numeric(summary(as.data.frame(table(network$TFSymbol))[, 
                                                                    2])[3])
  cut3 <- as.numeric(summary(as.data.frame(table(network$TFSymbol))[, 
                                                                    2])[5])
  vertex.size1 <- c()
  for (i in 1:nrow(nodes)) {
    if (nodes[i, 1] %in% tf_degree[, 1]) {
      if (tf_degree[tf_degree[, 1] == nodes[i, 1], ][, 
                                                     2] <= cut1) {
        vertex.size1 <- c(vertex.size1, vertex.size)
      }
      else if (tf_degree[tf_degree[, 1] == nodes[i, 1], 
      ][, 2] >= cut3) {
        vertex.size1 <- c(vertex.size1, vertex.size + 
                            (vertex.size.add) * 2)
      }
      else if (tf_degree[tf_degree[, 1] == nodes[i, 1], 
      ][, 2] < cut3 & tf_degree[tf_degree[, 1] == nodes[i, 
                                                        1], ][, 2] > cut1) {
        vertex.size1 <- c(vertex.size1, vertex.size + 
                            vertex.size.add)
      }
    }
    else {
      vertex.size1 <- c(vertex.size1, vertex.size)
    }
  }
  if (layout == "grid") {
    layout1 <- igraph::layout_on_grid(g)
  }
  else if (layout == "sphere") {
    layout1 <- igraph::layout_on_sphere(g)
  }
  else if (layout == "circle") {
    layout1 <- igraph::layout_in_circle(g)
  }
  else if (layout == "random") {
    layout1 <- igraph::layout_randomly(g)
  }
  else {
    stop("please input correct layout name")
  }
  edge.color2 <- c()
  for (i in edge_type) {
    if (i == "Positive") {
      edge.color2 <- c(edge.color2, edge.color[1])
    }
    else if (i == "Negative") {
      edge.color2 <- c(edge.color2, edge.color[2])
    }
  }
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$label.color <- edge.label.color
  igraph::E(g)$color <- edge.color2
  igraph::E(g)$width <- edge.width
  igraph::V(g)$color <- col1[nodes$type]
  igraph::V(g)$size <- vertex.size1
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$frame.color <- frame.color
  plot(g, layout = layout1, edge.curved = edge.curved, vertex.label.cex = vertex.label.cex, 
       layout = layout1, vertex.shape = "circle", vertex.label.family = vertex.label.family)
  if (legend == TRUE) {
    legend(x = 1.3, y = 1.3, paste0("Module", levels(factor(igraph::V(g)$type))), 
           pch = 21, col = "#777777", pt.bg = col1)
  }
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
}