library(WGCNA)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)
library(limma)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()
dist <- stats::dist

work_dir <- getwd()
setwd(work_dir)

expFilePath <- "merged_after_batch_removal.csv"
clusterFilePath <- "geneCluster.txt"
min_module_size <- 50
merge_threshold <- 0.25
soft_threshold_range <- 1:30
figure_width <- 12
figure_height <- 8

clusterData <- read.table(clusterFilePath, sep = "\t", header = TRUE, check.names = FALSE)
colnames(clusterData) <- c("ID", "ClusterType")
unique_clusters <- sort(unique(clusterData$ClusterType))

tmp_head <- readLines(expFilePath, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
rawData <- read.table(expFilePath, sep = sep, header = TRUE, check.names = FALSE)
rownames(rawData) <- rawData[, 1]
exprData <- rawData[, -1]
numericData <- matrix(as.numeric(as.matrix(exprData)), nrow = nrow(exprData),
                      dimnames = list(rownames(exprData), colnames(exprData)))
if (max(numericData) < 50) numericData <- log2(numericData + 1)
normalizedData <- normalizeBetweenArrays(numericData)
gene_sd <- apply(normalizedData, 1, sd, na.rm = TRUE)
filteredData <- normalizedData[gene_sd > 0.25, ]

matched_samples <- intersect(colnames(filteredData), clusterData$ID)
filteredData <- filteredData[, matched_samples]
clusterData <- clusterData[clusterData$ID %in% matched_samples, ]
clusterData <- clusterData[match(colnames(filteredData), clusterData$ID), ]

clinicalData <- data.frame(C1 = as.numeric(clusterData$ClusterType == "C1"),
                           C2 = as.numeric(clusterData$ClusterType == "C2"))
rownames(clinicalData) <- clusterData$ID

exprMatrix <- t(filteredData)

sftResult <- pickSoftThreshold(exprMatrix, powerVector = soft_threshold_range, RsquaredCut = 0.9, verbose = 0)
optimalPower <- sftResult$powerEstimate
if (is.na(optimalPower)) optimalPower <- 6

adjacencyMatrix <- adjacency(exprMatrix, power = optimalPower)
TOMMatrix <- TOMsimilarity(adjacencyMatrix)
dissTOM <- 1 - TOMMatrix
geneDendro <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(dendro = geneDendro, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = min_module_size)
moduleColors <- labels2colors(dynamicMods)
MEs0 <- moduleEigengenes(exprMatrix, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

merge <- mergeCloseModules(exprMatrix, moduleColors, cutHeight = merge_threshold, verbose = 0)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
moduleColors <- mergedColors
MEs <- mergedMEs

pdf("Figure8A.pdf", width = figure_width + 5, height = figure_height)
originalColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneDendro, cbind(originalColors, mergedColors),
                    c("Original", "Merged"), dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Module Assignment: Before and After Merging")
dev.off()

moduleTraitCor <- cor(MEs, clinicalData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))
moduleTraitPvalue_adj <- matrix(p.adjust(as.vector(moduleTraitPvalue), method = "BH"),
                                nrow = nrow(moduleTraitPvalue),
                                ncol = ncol(moduleTraitPvalue),
                                dimnames = dimnames(moduleTraitPvalue))
write.csv(data.frame(Module = rownames(moduleTraitCor), moduleTraitCor, check.names = FALSE),
          file = "Module_Trait_Correlation.csv", row.names = FALSE)
write.csv(data.frame(Module = rownames(moduleTraitPvalue_adj), moduleTraitPvalue_adj, check.names = FALSE),
          file = "Module_Trait_BHadjustedP.csv", row.names = FALSE)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue_adj, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf("Figure8B.pdf", width = 8, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinicalData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               cex.lab.x = 1.2,
               cex.lab.y = 1.0,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()
geneModules <- data.frame(Gene = colnames(exprMatrix), Module = moduleColors)
write.csv(geneModules, file = "Gene_Module_Assignment.csv", row.names = FALSE)

dotplot_enrich <- function(res, title, out_file, n_show = 15) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(invisible(NULL))
  df <- as.data.frame(res)
  df <- df[order(df$p.adjust), ]
  df <- df[seq_len(min(n_show, nrow(df))), ]
  p <- ggplot(df, aes(x = Count, y = reorder(Description, Count),
                      color = p.adjust, size = Count)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = title, x = "Gene Count", y = "", color = "p.adjust", size = "Count") +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(out_file, p, width = 10, height = 8)
}

panel_letters <- list(
  C1 = c(GO_BP = "C", GO_CC = "D", GO_MF = "E", KEGG = "F"),
  C2 = c(GO_BP = "G", GO_CC = "H", GO_MF = "I", KEGG = "J")
)

run_enrich <- function(genes, tag) {
  letters <- panel_letters[[tag]]
  if (is.null(letters)) letters <- c(GO_BP = "C", GO_CC = "D", GO_MF = "E", KEGG = "F")
  entrez <- bitr(unique(genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  if (length(entrez) == 0) return(invisible(NULL))
  go_bp <- tryCatch(enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
                    error = function(e) NULL)
  go_cc <- tryCatch(enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "CC",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
                    error = function(e) NULL)
  go_mf <- tryCatch(enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "MF",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
                    error = function(e) NULL)
  kegg <- tryCatch(enrichKEGG(gene = entrez, organism = "hsa",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                   error = function(e) NULL)
  if (!is.null(kegg)) kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  dotplot_enrich(go_bp, paste0(tag, " GO-BP"), paste0("Figure8", letters[["GO_BP"]], ".pdf"))
  dotplot_enrich(go_cc, paste0(tag, " GO-CC"), paste0("Figure8", letters[["GO_CC"]], ".pdf"))
  dotplot_enrich(go_mf, paste0(tag, " GO-MF"), paste0("Figure8", letters[["GO_MF"]], ".pdf"))
  dotplot_enrich(kegg, paste0(tag, " KEGG"), paste0("Figure8", letters[["KEGG"]], ".pdf"))
  if (!is.null(go_bp)) write.csv(as.data.frame(go_bp), paste0(tag, "_GO_BP.csv"), row.names = FALSE)
  if (!is.null(go_cc)) write.csv(as.data.frame(go_cc), paste0(tag, "_GO_CC.csv"), row.names = FALSE)
  if (!is.null(go_mf)) write.csv(as.data.frame(go_mf), paste0(tag, "_GO_MF.csv"), row.names = FALSE)
  if (!is.null(kegg)) write.csv(as.data.frame(kegg), paste0(tag, "_KEGG.csv"), row.names = FALSE)
}

for (cl in unique_clusters) {
  trait_col <- cl
  if (!(trait_col %in% colnames(moduleTraitCor))) next
  best_module <- rownames(moduleTraitCor)[which.max(moduleTraitCor[, trait_col])]
  best_color <- sub("^ME", "", best_module)
  module_genes <- geneModules$Gene[geneModules$Module == best_color]
  writeLines(module_genes, paste0(cl, "_module_genes.txt"))
  run_enrich(module_genes, cl)
}