library(limma)
library(ConsensusClusterPlus)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(reshape2)
library(ggpubr)
library(dplyr)
library(broom)
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(MCPcounter)

work_dir <- getwd()
setwd(work_dir)

gene_exp_file <- "merged_after_batch_removal.csv"
gene_list_file <- file.path(work_dir, "data", "hub_genes_9.csv")
cluster_file <- "geneCluster.txt"
cytokine_file <- file.path(work_dir, "data", "Inflammatory Cytokines gene.txt")
cli_file <- file.path(work_dir, "data", "clinical.xlsx")
disease_suffix <- "_tre"
maximum_k <- 11
final_cluster_num <- 2
total_iterations <- 50
item_sampling <- 0.8
feature_sampling <- 1
clust_algorithm <- "km"
distance_metric <- "euclidean"
consensus_seed <- 2026
set_rngvalue <- 2026
labelmap <- c("C1", "C2", "C3", "C4", "C5", "C6")

set.seed(set_rngvalue)
raw_data <- read.table(gene_exp_file, header = TRUE, sep = ",", check.names = FALSE)
gene_list <- read.csv(gene_list_file, header = TRUE, stringsAsFactors = FALSE)
target_genes <- gene_list[, 1]
raw_matrix <- as.matrix(raw_data)
rownames(raw_matrix) <- as.character(raw_matrix[, 1])
expr_data <- raw_matrix[, -1, drop = FALSE]
tre_cols <- grep(paste0(disease_suffix, "$"), colnames(expr_data), value = TRUE)
expr_data <- expr_data[, tre_cols, drop = FALSE]
matched_genes <- intersect(target_genes, rownames(expr_data))
expr_data <- expr_data[matched_genes, , drop = FALSE]
temp_matrix <- suppressWarnings(matrix(as.numeric(expr_data), nrow = nrow(expr_data)))
dimnames(temp_matrix) <- list(rownames(expr_data), colnames(expr_data))
final_data <- avereps(temp_matrix)

cluster_result <- ConsensusClusterPlus(final_data,
                                       maxK = maximum_k,
                                       reps = total_iterations,
                                       pItem = item_sampling,
                                       pFeature = feature_sampling,
                                       title = file.path(work_dir, "ConsensusClusterPlus"),
                                       clusterAlg = clust_algorithm,
                                       distance = distance_metric,
                                       seed = consensus_seed,
                                       plot = "png")

cls_vector <- cluster_result[[final_cluster_num]][["consensusClass"]]
cls_df <- data.frame(ClusterType = cls_vector)
uniq_levels <- levels(factor(cls_df$ClusterType))
cls_df$ClusterType <- labelmap[match(cls_df$ClusterType, uniq_levels)]
out_df <- rbind(ID = colnames(cls_df), cls_df)
write.table(out_df, file = cluster_file, sep = "\t", quote = FALSE, col.names = FALSE)
write.table(data.frame(Sample = rownames(cls_df), Cluster = cls_df$ClusterType),
            file = "Cluster_Labels.txt", sep = "\t", quote = FALSE, row.names = FALSE)

idx <- seq(0, 1, length.out = 101)
cdf_list <- list()
area_vec <- c()
for (k in 2:maximum_k) {
  cm <- cluster_result[[k]]$consensusMatrix
  fn <- ecdf(cm[upper.tri(cm)])
  cdf_list[[as.character(k)]] <- data.frame(k = k, x = idx, y = fn(idx))
  area_vec[as.character(k)] <- sum(diff(idx) * (head(fn(idx), -1) + tail(fn(idx), -1)) / 2)
}
cdf_df <- do.call(rbind, cdf_list)
cdf_df$k <- factor(cdf_df$k)
pdf("Figure7A.pdf", width = 7, height = 6)
print(ggplot(cdf_df, aes(x = x, y = y, color = k)) +
        geom_line(linewidth = 0.9) +
        labs(x = "Consensus index", y = "CDF", color = "k", title = "CDF curves") +
        theme_bw(base_size = 13))
dev.off()

pdf("Figure7B.pdf", width = 7, height = 6)
pheatmap(cluster_result[[final_cluster_num]]$consensusMatrix,
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("white", "#2166AC"))(100),
         main = paste0("Consensus matrix (k=", final_cluster_num, ")"))
dev.off()

track_mat <- sapply(2:maximum_k, function(k) cluster_result[[k]]$consensusClass)
track_mat <- track_mat[order(track_mat[, ncol(track_mat)]), , drop = FALSE]
track_df <- data.frame(Sample = factor(rownames(track_mat), levels = rownames(track_mat)),
                       k = factor(rep(2:maximum_k, each = nrow(track_mat))),
                       Cluster = factor(as.vector(track_mat)))
pdf("Figure7C.pdf", width = 9, height = 5)
print(ggplot(track_df, aes(x = Sample, y = k, fill = Cluster)) +
        geom_tile() +
        labs(x = "Sample", y = "k", title = "Tracking plot") +
        theme_minimal(base_size = 12) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid = element_blank()))
dev.off()

delta_df <- data.frame(k = 2:maximum_k, delta = c(NA, -diff(area_vec)))
pdf("Figure7D.pdf", width = 7, height = 6)
print(ggplot(delta_df, aes(x = k, y = delta)) +
        geom_line(linewidth = 0.9, color = "#2166AC") +
        geom_point(size = 3) +
        labs(x = "k", y = "Relative change in area under CDF", title = "Delta area") +
        theme_bw(base_size = 13))
dev.off()
expr <- final_data
group <- factor(cls_df[colnames(expr), "ClusterType"])
levs <- levels(group)
design <- model.matrix(~0 + group)
colnames(design) <- levs
fit <- lmFit(expr, design)
if (length(levs) == 2) {
  contrast_matrix <- makeContrasts(contrasts = paste0(levs[2], " - ", levs[1]), levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  deg_table <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "P")
} else {
  pair_contrasts <- combn(levs, 2, function(x) paste0(x[2], " - ", x[1]))
  contrast_matrix <- makeContrasts(contrasts = pair_contrasts, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  deg_table <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "F")
}
write.table(deg_table, file = "Cluster_DEG_results.txt", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = NA)

hub_expr <- data.frame(t(expr[matched_genes, , drop = FALSE]))
hub_expr$Cluster <- group
hub_long <- melt(hub_expr, id.vars = "Cluster",
                 variable.name = "Gene", value.name = "Expression")
hub_cor <- matrix(c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62"), ncol = 1)[1:length(levs), 1]
names(hub_cor) <- levs
p_hub <- ggplot(hub_long, aes(x = Gene, y = Expression, fill = Cluster)) +
  geom_boxplot(outlier.size = 1, width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = hub_cor) +
  labs(x = "", y = "Expression", fill = "Cluster") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave("Figure7E.pdf", p_hub, width = max(7, 0.6 * length(matched_genes) + 3), height = 5.5, dpi = 300)

cli <- as.data.frame(read_excel(cli_file))
rownames(cli) <- cli[[1]]
cli <- cli[, -1, drop = FALSE]
cluster_ids <- rownames(cls_df)
cluster_base <- gsub("_tre$", "", cluster_ids)
common_base <- intersect(cluster_base, rownames(cli))
keep_idx <- which(cluster_base %in% common_base)
cls_keep <- cls_df[keep_idx, , drop = FALSE]
cluster_base <- cluster_base[keep_idx]
cli <- cli[cluster_base, , drop = FALSE]
rownames(cli) <- rownames(cls_keep)

age_median <- median(as.numeric(cli$Age), na.rm = TRUE)
cli$Age <- ifelse(as.numeric(cli$Age) <= age_median, paste0("<=", age_median), paste0(">", age_median))
Type <- cbind(cls_keep, cli)
Type$ClusterType <- paste0("C", gsub("C", "", Type$ClusterType))
p_to_star <- function(p) {
  ifelse(is.na(p), "",
    ifelse(p < 0.001, "***",
      ifelse(p < 0.01, "**",
        ifelse(p < 0.05, "*", ""))))
}

clinical_vars <- colnames(Type)[-1]
chi_res <- data.frame()
for (clinical in clinical_vars) {
  tmp <- Type[c("ClusterType", clinical)]
  colnames(tmp) <- c("Cluster", "clinical")
  tmp <- tmp[!is.na(tmp[, "clinical"]) & tmp[, "clinical"] != "unknow", , drop = FALSE]
  st <- tryCatch(chisq.test(table(tmp)), error = function(e) list(p.value = NA_real_))
  chi_res <- rbind(chi_res, data.frame(Variable = clinical, P.Value = st$p.value))
}
chi_res$P.adj <- p.adjust(chi_res$P.Value, method = "BH", n = sum(!is.na(chi_res$P.Value)))
chi_res$Significance <- format_pval <- p_to_star(chi_res$P.adj)
write.csv(chi_res, file = "Cluster_Clinical_ChiSquare_Statistics.csv", row.names = FALSE)

sigVec <- c("ClusterType")
for (i in seq_len(nrow(chi_res))) {
  sigVec <- c(sigVec, paste0(chi_res$Variable[i], chi_res$Significance[i]))
}
colnames(Type) <- sigVec

exprRaw2 <- read.table(gene_exp_file, header = TRUE, sep = ",", check.names = FALSE)
exprMat2 <- as.matrix(exprRaw2)
rownames(exprMat2) <- exprMat2[, 1]
exprVals2 <- exprMat2[, -1, drop = FALSE]
numExpr2 <- matrix(as.numeric(as.matrix(exprVals2)), nrow = nrow(exprVals2),
                   dimnames = list(rownames(exprVals2), colnames(exprVals2)))
numExpr2 <- avereps(numExpr2)
data_hm <- numExpr2[matched_genes, rownames(Type), drop = FALSE]

colorList <- list()
nColors_total <- sum(sapply(colnames(Type), function(x) length(levels(factor(Type[, x])))))
bioCol <- colorRampPalette(brewer.pal(12, "Set3"))(nColors_total)
j <- 0
for (col_name in colnames(Type)) {
  cliLength <- length(levels(factor(Type[, col_name])))
  cliCol <- bioCol[(j + 1):(j + cliLength)]
  j <- j + cliLength
  names(cliCol) <- levels(factor(Type[, col_name]))
  if ("unknow" %in% names(cliCol)) cliCol["unknow"] <- "grey75"
  colorList[[col_name]] <- cliCol
}

pdf("Figure7F.pdf", width = 10, height = 6)
pheatmap(data_hm,
         annotation_col = Type,
         annotation_colors = colorList,
         color = colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "white",
                                    "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(100),
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         fontsize = 8,
         fontsize_row = 10,
         fontsize_col = 6,
         border_color = NA)
dev.off()

cli_raw <- as.data.frame(read_excel(cli_file))
rownames(cli_raw) <- cli_raw[[1]]
cli_raw <- cli_raw[, -1, drop = FALSE]
cli_raw <- cli_raw[cluster_base, , drop = FALSE]
rownames(cli_raw) <- rownames(cls_keep)
cli_raw$Cluster <- Type$ClusterType

all_var_names <- setdiff(colnames(cli_raw), "Cluster")
stats_box <- data.frame()
for (var_name in all_var_names) {
  vals <- cli_raw[[var_name]]
  is_numeric <- all(!is.na(suppressWarnings(as.numeric(vals))))
  box_vals <- if (is_numeric) as.numeric(vals) else as.numeric(factor(vals))
  box_tmp <- data.frame(Cluster = cli_raw$Cluster, Value = box_vals)
  box_tmp <- box_tmp[!is.na(box_tmp$Value), , drop = FALSE]
  n_groups <- length(unique(box_tmp$Cluster))
  if (n_groups < 2) next
  if (n_groups == 2) {
    wt <- tryCatch(wilcox.test(Value ~ Cluster, data = box_tmp), error = function(e) list(p.value = NA_real_))
    test_name <- "Wilcoxon rank-sum"
  } else {
    wt <- tryCatch(kruskal.test(Value ~ Cluster, data = box_tmp), error = function(e) list(p.value = NA_real_))
    test_name <- "Kruskal-Wallis"
  }
  stats_box <- rbind(stats_box, data.frame(Variable = var_name, Test = test_name, P.Value = wt$p.value))
}
stats_box$P.adj <- p.adjust(stats_box$P.Value, method = "BH", n = sum(!is.na(stats_box$P.Value)))
stats_box$Significance <- p_to_star(stats_box$P.adj)
write.csv(stats_box, file = "Cluster_Clinical_Boxplot_Statistics.csv", row.names = FALSE)

box_plot_list <- list()
for (var_name in all_var_names) {
  vals <- cli_raw[[var_name]]
  is_numeric <- all(!is.na(suppressWarnings(as.numeric(vals))))
  box_vals <- if (is_numeric) as.numeric(vals) else as.numeric(factor(vals))
  box_tmp <- data.frame(Cluster = cli_raw$Cluster, Value = box_vals)
  box_tmp <- box_tmp[!is.na(box_tmp$Value), , drop = FALSE]
  n_groups <- length(unique(box_tmp$Cluster))
  if (n_groups < 2) next
  padj <- stats_box$P.adj[match(var_name, stats_box$Variable)]
  label <- ifelse(is.na(padj), "NA",
            ifelse(padj < 0.001, "BH-adjusted P < 0.001",
                   paste0("BH-adjusted P = ", formatC(padj, format = "f", digits = 3))))
  box_palette <- c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62")[1:n_groups]
  p_box <- ggboxplot(box_tmp, x = "Cluster", y = "Value", fill = "Cluster",
                     palette = box_palette, xlab = "", ylab = var_name,
                     width = 0.5, notch = FALSE) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5, color = "red"),
          axis.text.x = element_text(size = 12, color = "black", face = "bold"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.line = element_line(color = "black", linewidth = 0.8),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10)) +
    labs(title = var_name, subtitle = label)
  box_plot_list[[var_name]] <- p_box
  ggsave(paste0("Figure7G_", var_name, ".pdf"), p_box, width = 4, height = 5)
}
if (length(box_plot_list) > 0) {
  combined_box <- ggarrange(plotlist = box_plot_list, ncol = length(box_plot_list), nrow = 1)
  ggsave("Figure7G.pdf", combined_box, width = 3.5 * length(box_plot_list), height = 5.5)
}
immune_abundance <- MCPcounter.estimate(
  expression = numExpr2,
  featuresType = "HUGO_symbols",
  genes = read.table(file.path(work_dir, "data", "genes.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
  probesets = read.table(file.path(work_dir, "data", "probesets.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
)
immune_df <- as.data.frame(t(immune_abundance))
write.csv(immune_df, file = "Figure7_MCPcounter_Results.csv", row.names = TRUE)

cluster_info <- read.table(cluster_file, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
colnames(cluster_info) <- "Cluster"
cluster_info$Cluster <- paste0("C", gsub("C", "", cluster_info$Cluster))
common_samples <- intersect(rownames(immune_df), rownames(cluster_info))
imm <- immune_df[common_samples, , drop = FALSE]
imm$Group <- cluster_info[common_samples, "Cluster"]
imm$Sample <- rownames(imm)
group_levels <- sort(unique(imm$Group))
n_groups <- length(group_levels)
group_colors <- c("#4393C3", "#D6604D", "#66C2A5", "#FC8D62")
names(group_colors) <- group_levels[1:min(n_groups, 4)]

imm_long <- melt(imm, id.vars = c("Sample", "Group"), variable.name = "Immune", value.name = "Abundance")
if (n_groups == 2) {
  imm_p <- imm_long %>% filter(!is.na(Group)) %>% group_by(Immune) %>%
    do(tidy(wilcox.test(Abundance ~ Group, data = .))) %>% dplyr::select(Immune, p.value)
} else {
  imm_p <- imm_long %>% filter(!is.na(Group)) %>% group_by(Immune) %>%
    do(tidy(kruskal.test(Abundance ~ Group, data = .))) %>% dplyr::select(Immune, p.value)
}
imm_p <- imm_p %>% ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         Sig = ifelse(p.adj < 0.001, "***", ifelse(p.adj < 0.01, "**", ifelse(p.adj < 0.05, "*", "ns"))))
write.csv(imm_p, file = "Immune_Cluster_PValues.csv", row.names = FALSE)

sig_y <- imm_long %>% filter(!is.na(Group)) %>% group_by(Immune) %>%
  summarise(y = max(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y + 0.05 * diff(range(imm_long$Abundance, na.rm = TRUE)))
sig_df <- imm_p %>% select(Immune, Sig) %>% left_join(sig_y, by = "Immune") %>%
  mutate(Immune = factor(Immune, levels = levels(imm_long$Immune)))

p_imm <- ggboxplot(imm_long, x = "Immune", y = "Abundance", fill = "Group",
                   palette = group_colors[1:n_groups],
                   xlab = "", ylab = "Abundance Score", legend.title = "Cluster",
                   notch = FALSE, width = 0.8) +
  geom_text(data = sig_df, aes(x = Immune, y = y, label = Sig),
            inherit.aes = FALSE, vjust = 0, size = 5, fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        plot.margin = margin(20, 10, 10, 10)) +
  coord_cartesian(clip = "off")
ggsave("Figure7H.pdf", p_imm, width = max(8, n_groups * 2.5), height = 6)

cytoGenes <- readLines(cytokine_file)
cytoGenes <- trimws(cytoGenes)
cytoGenes <- cytoGenes[cytoGenes != ""]
cytoMatched <- cytoGenes[cytoGenes %in% rownames(numExpr2)]
common_cyto <- intersect(colnames(numExpr2), rownames(cluster_info))
cyto_df <- as.data.frame(t(numExpr2[cytoMatched, common_cyto, drop = FALSE]))
cyto_df$Cluster <- cluster_info[common_cyto, "Cluster"]
cyto_df$Sample <- rownames(cyto_df)
cyto_long <- melt(cyto_df, id.vars = c("Sample", "Cluster"),
                  variable.name = "Cytokine", value.name = "Expression")

if (n_groups == 2) {
  cyto_p <- cyto_long %>% group_by(Cytokine) %>%
    do(tidy(wilcox.test(Expression ~ Cluster, data = .))) %>% dplyr::select(Cytokine, p.value)
} else {
  cyto_p <- cyto_long %>% group_by(Cytokine) %>%
    do(tidy(kruskal.test(Expression ~ Cluster, data = .))) %>% dplyr::select(Cytokine, p.value)
}
cyto_p <- cyto_p %>% ungroup() %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"),
         Sig = ifelse(p.adj < 0.001, "***", ifelse(p.adj < 0.01, "**", ifelse(p.adj < 0.05, "*", "ns"))))
write.csv(cyto_p, file = "Cytokine_Cluster_PValues.csv", row.names = FALSE)

sig_y_c <- cyto_long %>% group_by(Cytokine) %>%
  summarise(y = max(Expression, na.rm = TRUE), .groups = "drop") %>%
  mutate(y = y + 0.05 * diff(range(cyto_long$Expression, na.rm = TRUE)))
sig_c <- cyto_p %>% select(Cytokine, Sig) %>% left_join(sig_y_c, by = "Cytokine") %>%
  mutate(Cytokine = factor(Cytokine, levels = levels(cyto_long$Cytokine)))

p_cyto <- ggboxplot(cyto_long, x = "Cytokine", y = "Expression", fill = "Cluster",
                    palette = group_colors[1:n_groups],
                    xlab = "", ylab = "Expression Level", legend.title = "Cluster",
                    notch = FALSE, width = 0.8) +
  geom_text(data = sig_c, aes(x = Cytokine, y = y, label = Sig),
            inherit.aes = FALSE, vjust = 0, size = 5, fontface = "bold") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks = element_line(color = "black", linewidth = 1),
        plot.margin = margin(20, 10, 10, 10)) +
  coord_cartesian(clip = "off")
ggsave("Figure7I.pdf", p_cyto, width = max(10, 0.5 * length(cytoMatched)), height = 6)