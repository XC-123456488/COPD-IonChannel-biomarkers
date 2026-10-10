library(data.table)
library(limma)
library(readxl)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(broom)
library(ggsignif)
library(RColorBrewer)
library(patchwork)
library(tools)
library(ggrepel)
library(pheatmap)
library(grid)
library(ComplexHeatmap)
library(circlize)
library(ggvenn)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggprism)
library(ReactomePA)
library(reactome.db)

work_dir <- getwd()
setwd(work_dir)
raw_dir <- file.path(work_dir, "raw")
expr_dir <- work_dir

dataset_ids    <- c("GSE76925", "GSE20257", "GSE106986", "GSE38974", "GSE8545")
discovery_ids  <- c("GSE76925", "GSE20257")
gene_col_index <- 11

for (ds in dataset_ids) {
  ds_dir <- file.path(raw_dir, ds)
  if (!dir.exists(ds_dir)) next
  allFiles <- list.files(ds_dir, pattern = "\\.(txt|TXT)$")
  gseFiles <- allFiles[grepl("^GSE", allFiles, ignore.case = TRUE)]
  gplFiles <- allFiles[grepl("^GPL", allFiles, ignore.case = TRUE)]
  if (length(gseFiles) == 0 || length(gplFiles) == 0) next
  exprFilePath <- file.path(ds_dir, gseFiles[1])
  platformFilePath <- file.path(ds_dir, gplFiles[1])
  targetColIdx <- gene_col_index - 1
  allLines <- readLines(exprFilePath)
  idRefLine <- which(grepl("ID_REF", allLines, ignore.case = FALSE))[1]
  if (is.na(idRefLine)) stop(paste("no ID_REF line:", ds))
  exprData <- read.delim(exprFilePath, header = FALSE, sep = "\t", quote = "\"",
                         skip = idRefLine - 1, comment.char = "")
  colnames(exprData) <- exprData[1, ]
  exprData <- exprData[-1, ]
  colnames(exprData)[1] <- "ProbeID"
  rownames(exprData) <- NULL
  platformData <- read.delim(platformFilePath, header = FALSE, sep = "\t",
                             quote = "\"", comment.char = "#", stringsAsFactors = FALSE)
  geneMapping <- list()
  for (row in seq_len(nrow(platformData))) {
    if (all(platformData[row, ] == "") || grepl("^(ID|\\!)", platformData[row, 1])) next
    if (ncol(platformData) >= (targetColIdx + 1)) {
      rawGene <- platformData[row, targetColIdx + 1]
      if (rawGene != "" && !grepl(".+\\s+.+", rawGene)) {
        cleanGene <- sub("(.+?)///(.+)", "\\1", rawGene)
        cleanGene <- gsub('"', "", cleanGene)
        geneMapping[[as.character(platformData[row, 1])]] <- cleanGene
      }
    }
  }
  probeIDs <- exprData[["ProbeID"]]
  exprDataCore <- exprData[, -1, drop = FALSE]
  for (col in colnames(exprDataCore)) {
    exprDataCore[[col]] <- as.numeric(as.character(exprDataCore[[col]]))
  }
  exprDataWithID <- data.frame(ProbeID = probeIDs, exprDataCore, stringsAsFactors = FALSE)
  mappingDF <- data.frame(ProbeID = names(geneMapping),
                          geneSymbol = unlist(geneMapping),
                          stringsAsFactors = FALSE)
  mergedData <- merge(exprDataWithID, mappingDF, by = "ProbeID")
  sampleCols <- setdiff(colnames(mergedData), c("ProbeID", "geneSymbol"))
  geneGroups <- split(mergedData, mergedData$geneSymbol)
  geneNamesUnique <- names(geneGroups)
  resultList <- vector("list", length(geneNamesUnique))
  for (i in seq_along(geneNamesUnique)) {
    means <- colMeans(geneGroups[[i]][, sampleCols, drop = FALSE], na.rm = TRUE)
    resultList[[i]] <- c(geneSymbol = geneNamesUnique[i], means)
  }
  aggData <- as.data.frame(do.call(rbind, resultList), stringsAsFactors = FALSE)
  for (col in colnames(aggData)[-1]) aggData[[col]] <- as.numeric(aggData[[col]])
  aggData <- aggData[order(aggData$geneSymbol), ]
  write.csv(aggData, file = file.path(expr_dir, paste0(ds, "_geneMatrix.csv")), row.names = FALSE)
}
for (ds in dataset_ids) {
  f <- file.path(expr_dir, paste0(ds, "_geneMatrix.csv"))
  if (!file.exists(f)) next
  dt <- data.table::fread(f, header = TRUE, sep = ",")
  mat <- as.matrix(dt)
  rownames(mat) <- mat[, 1]
  expr <- mat[, -1, drop = FALSE]
  numericMat <- matrix(as.numeric(as.matrix(expr)), nrow = nrow(expr),
                       dimnames = list(rownames(expr), colnames(expr)))
  avgMat <- limma::avereps(numericMat)
  qx <- as.numeric(quantile(avgMat, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  needLog <- (qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0)
  if (needLog) {
    avgMat[avgMat < 0] <- 0
    avgMat <- log2(avgMat + 1)
  }
  normMat <- limma::normalizeBetweenArrays(avgMat)
  if (any(is.na(normMat))) normMat[is.na(normMat)] <- 0
  sampleInfo <- readxl::read_excel(file.path(raw_dir, ds, "sample_info.xlsx"), col_names = TRUE)
  if (ncol(sampleInfo) < 2) stop("sample_info.xlsx for ", ds, " must have two columns: sample ID and group label (con/tre)")
  sampleNames <- as.character(sampleInfo[[1]])
  sampleTypes <- as.character(sampleInfo[[2]])
  validIdx <- sampleNames %in% colnames(normMat)
  sampleNames <- sampleNames[validIdx]
  sampleTypes <- sampleTypes[validIdx]
  combinedData <- normMat[, sampleNames, drop = FALSE]
  sampleTypes <- ifelse(grepl("^(control|con)$", sampleTypes, ignore.case = TRUE), "con",
                    ifelse(grepl("^(disease|tre)$", sampleTypes, ignore.case = TRUE), "tre", sampleTypes))
  unknown_labels <- setdiff(unique(sampleTypes[!is.na(sampleTypes)]), c("con", "tre"))
  if (length(unknown_labels) > 0) stop("unrecognised group labels in sample_info.xlsx for ", ds, ": ", paste(unknown_labels, collapse = ", "), " (use con/tre)")
  colnames(combinedData) <- paste0(colnames(combinedData), "_", sampleTypes)
  finalMatrix <- cbind(GeneName = rownames(combinedData), combinedData)
  write.table(finalMatrix, file = file.path(expr_dir, paste0(ds, ".csv")),
              sep = ",", quote = FALSE, col.names = TRUE, row.names = FALSE)
}

hub_genes <- read.csv(file.path(work_dir, "data", "hub_genes_9.csv"), stringsAsFactors = FALSE)[[1]]
hub_genes <- trimws(hub_genes)
hub_genes <- hub_genes[hub_genes != ""]

for (ds in setdiff(dataset_ids, discovery_ids)) {
  f <- file.path(expr_dir, paste0(ds, ".csv"))
  if (!file.exists(f)) next
  dt <- data.table::fread(f)
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- dt[[1]]
  mode(mat) <- "numeric"
  keep <- hub_genes[hub_genes %in% rownames(mat)]
  if (length(keep) == 0) next
  out <- data.frame(Sample = colnames(mat), t(mat[keep, , drop = FALSE]), check.names = FALSE)
  write.csv(out, file = file.path(work_dir, "data", paste0("test_", ds, ".csv")), row.names = FALSE)
}

data_path <- expr_dir
file_list <- file.path(data_path, paste0(discovery_ids, ".csv"))
file_list <- file_list[file.exists(file_list)]
if (length(file_list) == 0) stop("no CSV files found for the discovery cohorts")
data_list <- lapply(file_list, function(f) {
  dt <- data.table::fread(f)
  data.table::setnames(dt, 1, "geneSymbol")
  dt
})
merged_data <- Reduce(function(x, y) merge(x, y, by = "geneSymbol"), data_list)
write.csv(merged_data, file = file.path(work_dir, "merged_before_batch_removal.csv"), row.names = FALSE)

batch_names <- tools::file_path_sans_ext(basename(file_list))
sample_counts <- sapply(data_list, function(x) ncol(x) - 1)
batch_info <- unlist(mapply(rep, batch_names, sample_counts))

expr_matrix <- as.data.frame(merged_data[, -1, with = FALSE])
rownames(expr_matrix) <- merged_data$geneSymbol
expr_matrix <- as.matrix(expr_matrix)
mode(expr_matrix) <- "numeric"
if (length(batch_info) != ncol(expr_matrix)) stop("batch vector length mismatch")

corrected_expr <- limma::removeBatchEffect(expr_matrix, batch = batch_info)
corrected_data <- data.frame(geneSymbol = rownames(expr_matrix), corrected_expr, check.names = FALSE)
write.csv(corrected_data, file = file.path(work_dir, "merged_after_batch_removal.csv"), row.names = FALSE)
expr_df <- as.data.frame(expr_matrix)
expr_df$gene <- rownames(expr_matrix)
expr_melt <- reshape2::melt(expr_df, id.vars = "gene", variable.name = "Sample", value.name = "Expression")
expr_melt$Batch <- batch_info[match(expr_melt$Sample, colnames(expr_matrix))]

n_batch <- length(unique(batch_info))
nature_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
if (n_batch > length(nature_colors)) nature_colors <- colorRampPalette(nature_colors)(n_batch)

p_before <- ggplot(expr_melt, aes(x = Sample, y = Expression, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, color = NA, lwd = 0) +
  scale_fill_manual(values = nature_colors[1:n_batch]) +
  labs(title = "Before Batch Effect Removal", x = "Sample", y = "Expression") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10))

corr_expr_df <- as.data.frame(corrected_expr)
corr_expr_df$gene <- rownames(expr_matrix)
corr_expr_melt <- reshape2::melt(corr_expr_df, id.vars = "gene", variable.name = "Sample", value.name = "Expression")
corr_expr_melt$Batch <- batch_info[match(corr_expr_melt$Sample, colnames(expr_matrix))]
p_after <- ggplot(corr_expr_melt, aes(x = Sample, y = Expression, fill = Batch)) +
  geom_boxplot(outlier.size = 0.3, color = NA, lwd = 0) +
  scale_fill_manual(values = nature_colors[1:n_batch]) +
  labs(title = "After Batch Effect Removal", x = "Sample", y = "Expression") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(color = "black", size = 0.8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.position = "right",
        plot.margin = margin(10, 10, 10, 10))

pca_plot <- function(mat, batch_info, title) {
  data_pca <- prcomp(t(mat), scale. = TRUE)
  pc_var <- data_pca$sdev^2 / sum(data_pca$sdev^2)
  pca_df <- data.frame(Sample = rownames(data_pca$x),
                       PC1 = data_pca$x[, 1],
                       PC2 = data_pca$x[, 2],
                       Batch = factor(batch_info, levels = unique(batch_info)))
  nb <- length(unique(batch_info))
  mycol <- brewer.pal(min(nb, 8), "Set1")
  if (nb > 8) mycol <- colorRampPalette(brewer.pal(9, "Set1"))(nb)
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 4, alpha = 0.90) +
    stat_ellipse(aes(fill = Batch), geom = "polygon", alpha = 0.14, linetype = 2, show.legend = FALSE) +
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    labs(title = title,
         x = paste0("PC1 (", round(pc_var[1] * 100, 1), "%)"),
         y = paste0("PC2 (", round(pc_var[2] * 100, 1), "%)"),
         color = "Batch") +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_line(colour = "gray90", linetype = 2),
          legend.title = element_text(face = "bold"),
          legend.background = element_rect(colour = "black", fill = NA, size = 0.25),
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(face = "bold"))
  p
}

pca_pre <- pca_plot(expr_matrix, batch_info, "PCA Before Batch Effect Removal")
pca_post <- pca_plot(corrected_expr, batch_info, "PCA After Batch Effect Removal")
combined_plot <- (p_before + p_after) / (pca_pre + pca_post) + plot_annotation(tag_levels = "A")
ggsave(file.path(work_dir, "Figure2A.pdf"), combined_plot, width = 16, height = 10)
raw_df <- read.table(file.path(work_dir, "merged_after_batch_removal.csv"), header = TRUE, sep = ",",
                     check.names = FALSE, stringsAsFactors = FALSE)
rownames(raw_df) <- raw_df[[1]]
expr_df2 <- raw_df[, -1, drop = FALSE]
expr_mat <- apply(expr_df2, 2, as.numeric)
rownames(expr_mat) <- rownames(raw_df)

geneset_dir <- file.path(work_dir, "data", "gene set")
txt_files <- list.files(geneset_dir, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)
if (length(txt_files) == 0) stop("no gene set files")

pathway_scores <- data.frame(row.names = colnames(expr_mat))
for (txt_file in txt_files) {
  pathway_name <- gsub("\\.txt$", "", basename(txt_file), ignore.case = TRUE)
  genes <- read.table(txt_file, header = TRUE, stringsAsFactors = FALSE)
  gene_vector <- as.character(genes[[1]])
  gene_vector <- gene_vector[gene_vector != "" & !is.na(gene_vector)]
  matched_genes <- intersect(gene_vector, rownames(expr_mat))
  if (length(matched_genes) >= 2) {
    pathway_scores[[pathway_name]] <- colMeans(expr_mat[matched_genes, ], na.rm = TRUE)
  } else if (length(matched_genes) == 1) {
    pathway_scores[[pathway_name]] <- expr_mat[matched_genes, ]
  }
}

pathway_scores$Group <- ifelse(
  grepl("_con$", rownames(pathway_scores), ignore.case = TRUE), "Control",
  ifelse(grepl("_tre$", rownames(pathway_scores), ignore.case = TRUE), "Disease", NA))
pathway_scores$Sample <- rownames(pathway_scores)
pathway_scores <- pathway_scores[!is.na(pathway_scores$Group), ]

pathway_cols <- setdiff(colnames(pathway_scores), c("Group", "Sample"))
data_long <- reshape2::melt(pathway_scores, id.vars = c("Sample", "Group"),
                            measure.vars = pathway_cols,
                            variable.name = "Pathway", value.name = "Score")
data_long$Score <- as.numeric(data_long$Score)
write.csv(pathway_scores, file = file.path(work_dir, "Pathway_Mean_Expression_Scores.csv"), row.names = FALSE)

p_to_star <- function(p, ns = "ns") {
  ifelse(is.na(p), ns,
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ns))))
}

pvalue_table <- data_long %>%
  group_by(Pathway) %>%
  do({
    groups <- unique(.$Group)
    if (length(groups) == 2) {
      tryCatch(tidy(wilcox.test(Score ~ Group, data = .)),
               error = function(e) data.frame(p.value = NA_real_))
    } else {
      data.frame(p.value = NA_real_)
    }
  }) %>%
  dplyr::select(Pathway, p.value) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH", n = sum(!is.na(p.value))),
         Significance = p_to_star(p.adj)) %>%
  arrange(p.value)
write.csv(pvalue_table, file = file.path(work_dir, "pathway_pvalues.csv"), row.names = FALSE)

countControl <- sum(pathway_scores$Group == "Control")
countDisease <- sum(pathway_scores$Group == "Disease")
natureColors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
names(natureColors) <- c("Control", "Disease")

for (pw in pathway_cols) {
  pw_data <- data_long[data_long$Pathway == pw, ]
  pw_adj_p <- pvalue_table$p.adj[match(pw, pvalue_table$Pathway)]
  pw_label <- p_to_star(pw_adj_p)
  boxplot_fig <- ggboxplot(pw_data, x = "Group", y = "Score", fill = "Group",
                           palette = natureColors, xlab = "", ylab = paste0(pw, " Score"),
                           legend.title = "Group", notch = FALSE, width = 0.6) +
    geom_signif(comparisons = list(c("Control", "Disease")),
                annotations = pw_label,
                textsize = 5, vjust = -0.2, tip_length = 0.02) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top",
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black", face = "bold"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, color = "black", face = "bold"),
          axis.line = element_line(color = "black", linewidth = 0.8),
          axis.ticks = element_line(color = "black", linewidth = 0.8),
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5, color = "black"),
          plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black"),
          plot.margin = margin(20, 10, 10, 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA)) +
    labs(title = paste0(pw, " Expression Comparison"),
         subtitle = paste0("Control n=", countControl, ", Disease n=", countDisease)) +
    coord_cartesian(clip = "off")
  ggsave(file.path(work_dir, paste0("Figure2B_", pw, ".pdf")), boxplot_fig, width = 3, height = 6, dpi = 300)
}

summary_table <- data_long %>%
  group_by(Pathway, Group) %>%
  summarise(MeanScore = mean(Score, na.rm = TRUE),
            MedianScore = median(Score, na.rm = TRUE),
            SD = sd(Score, na.rm = TRUE),
            Count = n(), .groups = "drop")
write.csv(summary_table, file = file.path(work_dir, "pathway_summary_table.csv"), row.names = FALSE)

threshold_logFC <- 0
threshold_adjP <- 0.05

file_expr <- file.path(work_dir, "merged_after_batch_removal.csv")
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
expr_raw <- read.table(file_expr, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
rownames(expr_raw) <- expr_raw[, 1]
expr_mat2 <- expr_raw[, -1, drop = FALSE]
expr_mat2 <- as.matrix(expr_mat2)
expr_mat2 <- apply(expr_mat2, 2, as.numeric)
rownames(expr_mat2) <- rownames(expr_raw)
colnames(expr_mat2) <- colnames(expr_raw)[-1]
combined_expr <- expr_mat2

sample_names <- colnames(expr_mat2)
group_info <- ifelse(grepl("_con$", sample_names, ignore.case = TRUE), "Control",
                     ifelse(grepl("_tre$", sample_names, ignore.case = TRUE), "Disease", "Unknown"))
if (any(group_info == "Unknown")) stop("unknown sample group")
num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Disease")

group_labels <- factor(group_info, levels = c("Control", "Disease"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Disease")
fit <- lmFit(expr_mat2, design_mat)
contrast_mat <- makeContrasts(Disease - Control, levels = design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)
all_diff_results <- topTable(fit2, adjust.method = "fdr", number = Inf)
write.table(cbind(Gene = rownames(all_diff_results), all_diff_results),
            file = file.path(work_dir, "DE_results.csv"), sep = ",", quote = FALSE, row.names = FALSE)

significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]
output_DEGs <- cbind(Gene = rownames(significant_DEGs), significant_DEGs)
SE <- tryCatch(abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])), error = function(e) NA)
output_DEGs <- cbind(output_DEGs, SE = SE)
output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors = FALSE)
output_DEGs$Regulation <- ifelse(as.numeric(output_DEGs$logFC) > 0, "Up", "Down")
desired_order <- c("Gene", "logFC", "Regulation", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]
write.csv(output_DEGs, file = file.path(work_dir, "DE_significant_genes.csv"), row.names = FALSE)

gene_set_dir <- file.path(work_dir, "data", "gene set")
gene_set_files <- list.files(gene_set_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(gene_set_files) == 0) stop("no gene set txt files")
deg_genes <- rownames(significant_DEGs)
all_gene_sets <- list()
all_intersect_results <- list()
for (gs_file in gene_set_files) {
  gs_name <- tools::file_path_sans_ext(basename(gs_file))
  gs_data <- read.table(gs_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
  if ("Genes" %in% colnames(gs_data)) {
    gs_genes <- unique(gs_data$Genes)
  } else if ("Gene" %in% colnames(gs_data)) {
    gs_genes <- unique(gs_data$Gene)
  } else {
    gs_genes <- unique(gs_data[[1]])
  }
  gs_genes <- gs_genes[!is.na(gs_genes) & gs_genes != ""]
  gs_intersect <- intersect(gs_genes, deg_genes)
  all_gene_sets[[gs_name]] <- gs_genes
  all_intersect_results[[gs_name]] <- gs_intersect
}
research_genes <- all_gene_sets[[length(all_gene_sets)]]
intersect_genes <- all_intersect_results[[length(all_intersect_results)]]

volcano_data <- all_diff_results
volcano_data$Gene <- rownames(volcano_data)
intersect_deg_data <- volcano_data[volcano_data$Gene %in% intersect_genes, ]
if (nrow(intersect_deg_data) <= 40) {
  label_genes_gradient <- intersect_deg_data$Gene
} else {
  up_intersect <- intersect_deg_data[intersect_deg_data$logFC > 0, ]
  up_intersect <- up_intersect[order(-(-log10(up_intersect$adj.P.Val))), ]
  down_intersect <- intersect_deg_data[intersect_deg_data$logFC < 0, ]
  down_intersect <- down_intersect[order(-(-log10(down_intersect$adj.P.Val))), ]
  up_label <- head(up_intersect$Gene, min(20, nrow(up_intersect)))
  down_label <- head(down_intersect$Gene, min(20, nrow(down_intersect)))
  label_genes_gradient <- unique(c(up_label, down_label))
}
volcano_data$label_gene <- ifelse(volcano_data$Gene %in% label_genes_gradient, volcano_data$Gene, "")

x_limit <- max(abs(volcano_data$logFC), na.rm = TRUE) * 1.1
y_max <- max(-log10(volcano_data$adj.P.Val), na.rm = TRUE) * 1.05
y_values <- -log10(volcano_data$adj.P.Val)
y_values <- y_values[is.finite(y_values)]
y_range <- range(y_values, na.rm = TRUE)
size_breaks <- pretty(y_range, n = 5)
size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max(y_values)]

up_count <- sum(volcano_data$logFC > threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)
down_count <- sum(volcano_data$logFC < -threshold_logFC & volcano_data$adj.P.Val < threshold_adjP)

gradient_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha = 0.85) +
  geom_point(data = subset(volcano_data, label_gene != ""),
             aes(size = -log10(adj.P.Val)),
             shape = 21, fill = NA, color = "black", stroke = 1.2) +
  scale_color_gradientn(
    colors = c("#2B83BA", "#5AAE61", "#A6D96A", "#FFFFBF", "#FEE08B", "#FDAE61", "#D7191C"),
    values = scales::rescale(c(-6, -4, -2, 0, 2, 4, 6)),
    limits = c(-x_limit, x_limit)) +
  scale_size_continuous(range = c(1, 5), breaks = size_breaks) +
  geom_vline(xintercept = c(-threshold_logFC, threshold_logFC), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(threshold_adjP), linetype = "dashed", color = "grey40") +
  annotate("text", x = -x_limit * 0.98, y = y_max * 0.95,
           label = paste0("Down: ", down_count), hjust = 0, size = 4, color = "grey20") +
  annotate("text", x = x_limit * 0.98, y = y_max * 0.95,
           label = paste0("Up: ", up_count), hjust = 1, size = 4, color = "grey20") +
  annotate("text", x = x_limit * 0.98, y = -log10(threshold_adjP),
           label = paste0("p = ", threshold_adjP),
           hjust = 1, vjust = -0.5, size = 3.5, color = "grey30") +
  scale_x_continuous(limits = c(-x_limit, x_limit), expand = c(0.02, 0)) +
  scale_y_continuous(limits = c(0, y_max), expand = c(0.02, 0)) +
  labs(title = "Volcano Plot", x = "avg_log2FC", y = "-log10(p_val)") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black", size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1.2),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
        legend.key = element_rect(fill = "white"),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.spacing.y = unit(0.3, "cm"),
        plot.margin = margin(15, 15, 10, 10)) +
  guides(color = guide_colorbar(title = "avg_log2FC", title.position = "top", title.hjust = 0.5,
                                barwidth = 1.2, barheight = 10, frame.colour = "black",
                                frame.linewidth = 0.5, ticks.colour = "black",
                                ticks.linewidth = 0.5, order = 1),
         size = guide_legend(title = "-log10(p_val)", title.position = "top", title.hjust = 0.5,
                             override.aes = list(alpha = 1, color = "grey30"), order = 2))
ggsave(file.path(work_dir, "Figure2C.pdf"), gradient_volcano, width = 9, height = 7.5, dpi = 300)

diff_gene_list <- data.frame(gene = rownames(significant_DEGs))
write.table(diff_gene_list, file = file.path(work_dir, "DEG_geneList.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
for (gs_name in names(all_gene_sets)) {
  gs_genes <- all_gene_sets[[gs_name]]
  gs_intersect <- all_intersect_results[[gs_name]]
  gene_list_venn <- list()
  gene_list_venn[[gs_name]] <- gs_genes
  gene_list_venn[["DEG"]] <- deg_genes
  venn_plot <- ggvenn(gene_list_venn,
                      show_percentage = TRUE,
                      stroke_color = "white",
                      stroke_size = 0.5,
                      fill_color = c("#1E90FF", "#FF8C00"),
                      set_name_color = c("#1E90FF", "#FF8C00"),
                      set_name_size = 6,
                      text_size = 4.5)
  ggsave(file.path(work_dir, "Figure2D.pdf"), venn_plot, width = 8, height = 8)
  write.csv(data.frame(Gene = gs_intersect),
            file = file.path(work_dir, "intersection_genes_13.csv"), row.names = FALSE)
  write.table(data.frame(Gene = gs_intersect),
              file = file.path(work_dir, "intersection_genes_13.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

intersect_gene_file <- file.path(work_dir, "intersection_genes_13.txt")
intersect_gene_data <- read.table(intersect_gene_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
plot_genes <- intersect_gene_data[[1]]
plot_genes <- plot_genes[plot_genes %in% rownames(combined_expr)]

heatmap_mat <- combined_expr[plot_genes, , drop = FALSE]
heatmap_scaled <- t(scale(t(heatmap_mat)))
ha_col <- HeatmapAnnotation(
  Group = factor(c(rep("Control", num_ctrl), rep("Disease", num_treat))),
  col = list(Group = c("Control" = "#2166AC", "Disease" = "#B2182B")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11))
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
ht_intersect <- Heatmap(
  heatmap_scaled,
  name = "Z-score",
  col = col_fun,
  top_annotation = ha_col,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 10, fontface = "italic"),
  column_title = "Intersection Genes Heatmap",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 11, fontface = "bold"),
                              labels_gp = gpar(fontsize = 9)))
pdf(file.path(work_dir, "Figure2F.pdf"), width = 10, height = max(4, length(plot_genes) * 0.45))
draw(ht_intersect)
dev.off()

group_vec <- factor(c(rep("Control", num_ctrl), rep("Disease", num_treat)), levels = c("Control", "Disease"))
box_list <- list()
for (g in plot_genes) {
  box_list[[g]] <- data.frame(Gene = g,
                              Expression = as.numeric(combined_expr[g, ]),
                              Group = group_vec)
}
box_df <- do.call(rbind, box_list)
box_df$Gene <- factor(box_df$Gene, levels = plot_genes)

p_vals <- sapply(plot_genes, function(g) {
  if (g %in% rownames(all_diff_results)) all_diff_results[g, "adj.P.Val"] else NA
})
format_pval <- function(p) {
  if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else return("ns")
}
p_labels <- sapply(p_vals, format_pval)
y_positions <- sapply(plot_genes, function(g) max(as.numeric(combined_expr[g, ])) * 1.08)
anno_df <- data.frame(Gene = factor(plot_genes, levels = plot_genes),
                      label = p_labels,
                      y = y_positions)

box_plot <- ggplot(box_df, aes(x = Gene, y = Expression, fill = Group)) +
  geom_boxplot(outlier.size = 0.8, width = 0.65, alpha = 0.85, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("Control" = "#2166AC", "Disease" = "#B2182B")) +
  geom_text(data = anno_df, aes(x = Gene, y = y, label = label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  labs(x = "", y = "Expression", fill = "Group") +
  theme_bw(base_size = 13) +
  theme(axis.text.x = element_text(size = 11, color = "black", face = "bold.italic",
                                   angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", size = 13),
        legend.position = "top",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1))
ggsave(file.path(work_dir, "Figure2E.pdf"), box_plot,
       width = max(6, length(plot_genes) * 1.2 + 2) - 5, height = 5.5, dpi = 300)
gene_data <- read.table(file.path(work_dir, "DEG_geneList.txt"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE)
gene_list <- gene_data[[1]] %>% as.character() %>% trimws() %>% discard(~ is.na(.x) || .x == "") %>% unique()
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db) %>%
  distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
entrez_ids <- unique(gene_ids$ENTREZID)

go_bp <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 1, readable = TRUE)
go_cc <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_mf <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

kegg_obj <- enrichKEGG(gene = entrez_ids, organism = "hsa", pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 1)
if (!is.null(kegg_obj)) kegg_obj <- setReadable(kegg_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

reactome_obj <- ReactomePA::enrichPathway(gene = entrez_ids, organism = "human",
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(reactome_obj)) reactome_obj <- setReadable(reactome_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

do_obj <- enrichDO(gene = entrez_ids, pAdjustMethod = "BH", pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2, readable = TRUE)
if (!is.null(do_obj)) do_obj <- setReadable(do_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

immport_file <- list.files(work_dir, pattern = "(?i)(immport|all_gene_lists).*[.](gmt|csv|tsv|txt|xlsx)$",
                           full.names = TRUE)
immport_file <- immport_file[!grepl("(?i)_enrich|result", immport_file, perl = TRUE)]
if (length(immport_file) > 0) {
  term2gene <- read.gmt(immport_file[1]) %>%
    transmute(term = as.character(term), gene = as.character(gene)) %>%
    separate_rows(gene, sep = "[,;/| ]+") %>%
    mutate(term = trimws(term), gene = trimws(gene)) %>%
    filter(!is.na(term), !is.na(gene), term != "", gene != "") %>%
    distinct(term, gene)
  immport_enrich <- enricher(gene = gene_list, TERM2GENE = term2gene, pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                             minGSSize = 5, maxGSSize = 5000)
} else {
  immport_enrich <- NULL
}

as_result_df <- function(enrich_obj, ontology = NULL) {
  if (is.null(enrich_obj)) return(data.frame())
  result_df <- as.data.frame(enrich_obj)
  if (nrow(result_df) == 0) return(data.frame())
  if (!is.null(ontology)) result_df <- result_df %>% mutate(ONTOLOGY = ontology)
  result_df
}

go_all <- bind_rows(as_result_df(go_bp, "BP"), as_result_df(go_cc, "CC"), as_result_df(go_mf, "MF"))
kegg_df <- as_result_df(kegg_obj, "KEGG")
reactome_df <- as_result_df(reactome_obj, "Reactome")
do_df <- as_result_df(do_obj, "DO")
immport_df <- as_result_df(immport_enrich, "ImmPort")

write.csv(go_all, file = file.path(work_dir, "GO_enrich.csv"), row.names = FALSE)
write.csv(kegg_df, file = file.path(work_dir, "KEGG_enrich.csv"), row.names = FALSE)
write.csv(reactome_df, file = file.path(work_dir, "Reactome_enrich.csv"), row.names = FALSE)
write.csv(do_df, file = file.path(work_dir, "DO_enrich.csv"), row.names = FALSE)
write.csv(immport_df, file = file.path(work_dir, "ImmPort_enrich.csv"), row.names = FALSE)

selected_list <- list()
if (nrow(go_all) > 0) {
  selected_list <- append(selected_list, list(go_all %>% group_by(ONTOLOGY) %>%
                                                arrange(p.adjust) %>% slice_head(n = 5) %>% ungroup()))
}
if (nrow(kegg_df) > 0) {
  selected_list <- append(selected_list, list(kegg_df %>% arrange(p.adjust) %>%
                                                slice_head(n = 5) %>% mutate(ONTOLOGY = "KEGG")))
}
if (nrow(reactome_df) > 0) {
  selected_list <- append(selected_list, list(reactome_df %>% arrange(p.adjust) %>%
                                                slice_head(n = 5) %>% mutate(ONTOLOGY = "Reactome")))
}
if (nrow(do_df) > 0) {
  selected_list <- append(selected_list, list(do_df %>% arrange(p.adjust) %>%
                                                slice_head(n = 5) %>% mutate(ONTOLOGY = "DO")))
}
if (nrow(immport_df) > 0) {
  selected_list <- append(selected_list, list(immport_df %>% arrange(p.adjust) %>%
                                                slice_head(n = 5) %>% mutate(ONTOLOGY = "ImmPort")))
}

ontology_level <- c("BP", "CC", "MF", "KEGG", "Reactome", "DO", "ImmPort")
pal_all <- c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087", KEGG = "#3C5488",
             Reactome = "#F39C12", DO = "#BE04B9", ImmPort = "#7E6148")

use_pathway <- bind_rows(selected_list) %>%
  filter(!is.na(p.adjust), p.adjust > 0) %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(ontology_level))) %>%
  arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = unique(Description)),
         geneID = map_chr(str_split(geneID, "/"), ~ paste(head(.x, 10), collapse = "/"))) %>%
  tibble::rowid_to_column("index")

width <- 1
xaxis_max <- max(-log10(use_pathway$p.adjust), na.rm = TRUE) + 5
xaxis_break_by <- ifelse(xaxis_max <= 10, 2, 10)
active_levels <- levels(droplevels(use_pathway$ONTOLOGY))
pal <- pal_all[rev(active_levels)]

rect.data <- use_pathway %>%
  count(ONTOLOGY, name = "n") %>%
  mutate(xmin = -4 * width,
         xmax = -2 * width,
         ymax = cumsum(n),
         ymin = lag(ymax, default = 0) + 0.6,
         ymax = ymax + 0.4)

p <- ggplot(use_pathway, aes(x = -log10(p.adjust), y = index, fill = ONTOLOGY)) +
  geom_col(width = 0.6, alpha = 0.8, orientation = "y") +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 5) +
  geom_text(aes(x = 0.1, label = geneID, colour = ONTOLOGY),
            hjust = 0, vjust = 2.4, size = 3.5, fontface = "italic", show.legend = FALSE) +
  geom_point(aes(x = -width, size = Count), shape = 21) +
  geom_text(aes(x = -width, label = Count), size = 3) +
  scale_size_continuous(name = "Count", range = c(5, 12)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
            data = rect.data, inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect.data, angle = 90, inherit.aes = FALSE) +
  annotate("segment", x = 0, xend = xaxis_max, y = 0, yend = 0, linewidth = 1.5) +
  labs(y = NULL, x = "-log10(p.adjust)") +
  scale_fill_manual(values = pal_all) +
  scale_colour_manual(values = pal_all) +
  scale_x_continuous(breaks = seq(0, xaxis_max, xaxis_break_by),
                     expand = expansion(mult = c(0.05, 0))) +
  guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE)) +
  theme_prism() +
  theme(axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks.y = element_blank())
ggsave(file.path(work_dir, "Figure2H.pdf"), p, width = 12, height = 13)