library(e1071)
library(parallel)
library(preprocessCore)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(broom)
library(ggsci)
library(limma)
library(pheatmap)
library(scales)
library(MCPcounter)

work_dir <- getwd()
setwd(work_dir)

input_file     <- "merged_after_batch_removal.csv"
refer_file     <- file.path(work_dir, "data", "refer.txt")
gene_list_file <- file.path(work_dir, "data", "hub_genes_9.csv")
cytokine_file  <- file.path(work_dir, "data", "Inflammatory Cytokines gene.txt")

bubble_plot <- function(df, y_col, ylab_title, out_file, w, h) {
  df$AbsCorr <- abs(df$Correlation)
  df$Sig <- ifelse(df$P.Adj < 0.001, "***",
            ifelse(df$P.Adj < 0.01, "**",
            ifelse(df$P.Adj < 0.05, "*", "")))
  df[[y_col]] <- factor(df[[y_col]], levels = rev(sort(unique(df[[y_col]]))))
  df$Gene <- factor(df$Gene, levels = sort(unique(as.character(df$Gene))))
  p <- ggplot(df, aes(x = Gene, y = .data[[y_col]])) +
    geom_point(aes(size = AbsCorr, fill = Correlation),
               shape = 22, color = "grey30", stroke = 0.3) +
    scale_fill_viridis_c(option = "viridis",
                         limits = c(-max(df$AbsCorr), max(df$AbsCorr)),
                         oob = squish, name = "Correlation",
                         breaks = pretty(range(df$Correlation), n = 5)) +
    scale_size_continuous(range = c(2, 10), guide = "none") +
    geom_text(data = df[df$Sig != "", ],
              aes(label = Sig), color = "white", size = 4, fontface = "bold", vjust = 0.5) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", size = 1.2),
      panel.grid.major = element_line(color = "black", size = 0.6),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1, face = "italic"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.8),
      legend.position = "right",
      legend.key.height = unit(1.8, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(10, 10, 10, 10)
    )
  pdf(out_file, width = w, height = h)
  print(p)
  dev.off()
}

coreAlgorithm <- function(refer_matrix, target_vector) {
  num_steps <- 3
  runSVM <- function(step_index) {
    nu_param <- c(0.25, 0.5, 0.75)[step_index]
    svm(refer_matrix, target_vector, type = "nu-regression", kernel = "linear", nu = nu_param, scale = FALSE)
  }
  if (Sys.info()["sysname"] == "Windows") {
    models_list <- mclapply(1:num_steps, runSVM, mc.cores = 1)
  } else {
    models_list <- mclapply(1:num_steps, runSVM, mc.cores = num_steps)
  }
  rmse_vec <- rep(0, num_steps)
  corr_vec <- rep(0, num_steps)
  for (i in 1:num_steps) {
    current_model <- models_list[[i]]
    weight_vals <- t(current_model$coefs) %*% current_model$SV
    weight_vals[which(weight_vals < 0)] <- 0
    norm_weights <- weight_vals / sum(weight_vals)
    weighted_refer <- sweep(refer_matrix, MARGIN = 2, norm_weights, "*")
    est_target <- apply(weighted_refer, 1, sum)
    rmse_vec[i] <- sqrt(mean((est_target - target_vector)^2))
    corr_vec[i] <- cor(est_target, target_vector)
  }
  best_idx <- which.min(rmse_vec)
  best_model <- models_list[[best_idx]]
  best_weights <- t(best_model$coefs) %*% best_model$SV
  best_weights[which(best_weights < 0)] <- 0
  final_weights <- best_weights / sum(best_weights)
  list(final_weights = final_weights, best_rmse = rmse_vec[best_idx], best_corr = corr_vec[best_idx])
}

doPermutation <- function(num_perm, refer_matrix, mix_matrix) {
  target_list <- as.list(data.matrix(mix_matrix))
  corr_distribution <- numeric()
  for (perm_counter in 1:num_perm) {
    permuted_target <- as.numeric(target_list[sample(length(target_list), nrow(refer_matrix), replace = TRUE)])
    if (sd(permuted_target) == 0) {
      std_target <- permuted_target
    } else {
      std_target <- (permuted_target - mean(permuted_target)) / sd(permuted_target)
    }
    perm_result <- coreAlgorithm(refer_matrix, std_target)
    corr_distribution <- c(corr_distribution, perm_result$best_corr)
  }
  list(null_distribution = sort(corr_distribution))
}

CIBERSORTModified <- function(refer_file, mix_file, num_perm = 0, do_QN = TRUE) {
  refer_data <- data.matrix(read.table(refer_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
  mix_data <- data.matrix(read.table(mix_file, header = TRUE, sep = ",", row.names = 1, check.names = FALSE))
  refer_data <- refer_data[order(rownames(refer_data)), ]
  mix_data <- mix_data[order(rownames(mix_data)), ]
  if (max(mix_data) < 50) mix_data <- 2^mix_data
  if (do_QN) {
    orig_colnames <- colnames(mix_data)
    orig_rownames <- rownames(mix_data)
    mix_data <- normalize.quantiles(mix_data)
    colnames(mix_data) <- orig_colnames
    rownames(mix_data) <- orig_rownames
  }
  common_genes <- intersect(rownames(refer_data), rownames(mix_data))
  refer_data <- refer_data[common_genes, , drop = FALSE]
  mix_data <- mix_data[common_genes, , drop = FALSE]
  refer_data <- (refer_data - mean(refer_data)) / sd(as.vector(refer_data))
  null_dist <- NULL
  if (num_perm > 0) null_dist <- doPermutation(num_perm, refer_data, mix_data)$null_distribution
  output_header <- c("Mixture", colnames(refer_data), "P-value", "Correlation", "RMSE")
  output_matrix <- NULL
  total_samples <- ncol(mix_data)
  for (sample_idx in 1:total_samples) {
    curr_sample <- mix_data[, sample_idx]
    std_sample <- if (sd(curr_sample) == 0) curr_sample else (curr_sample - mean(curr_sample)) / sd(curr_sample)
    alg_result <- coreAlgorithm(refer_data, std_sample)
    weights <- alg_result$final_weights
    if (!is.null(null_dist)) {
      p_value <- mean(null_dist >= alg_result$best_corr)
    } else {
      p_value <- 1
    }
    output_matrix <- rbind(output_matrix, c(colnames(mix_data)[sample_idx], weights, p_value, alg_result$best_corr, alg_result$best_rmse))
  }
  colnames(output_matrix) <- output_header
  rownames(output_matrix) <- colnames(mix_data)
  as.data.frame(output_matrix, stringsAsFactors = FALSE)
}
result_table <- CIBERSORTModified(refer_file, input_file, num_perm = 1000, do_QN = TRUE)
filtered_table <- result_table[result_table[, "P-value"] < 1, ]
weight_only <- as.matrix(filtered_table[, 1:(ncol(filtered_table) - 3)])
final_output <- rbind(id = colnames(weight_only), weight_only)
write.table(final_output, file = "CIBERSORT-Results.csv", sep = ",", quote = FALSE, col.names = FALSE)

cibersort_rt <- as.data.frame(weight_only)
cibersort_rt[] <- lapply(cibersort_rt, function(x) as.numeric(as.character(x)))

exprRaw <- read.table(input_file, header = TRUE, sep = ",", check.names = FALSE)
exprMat <- as.matrix(exprRaw)
rownames(exprMat) <- exprMat[, 1]
exprVals <- exprMat[, -1, drop = FALSE]
numExpr <- matrix(as.numeric(as.matrix(exprVals)), nrow = nrow(exprVals),
                  dimnames = list(rownames(exprVals), colnames(exprVals)))
numExpr <- avereps(numExpr)

gene_list <- read.csv(gene_list_file, header = TRUE, stringsAsFactors = FALSE)
targetGene <- as.character(gene_list[, 1])
targetGene <- targetGene[targetGene != ""]
matched_genes <- targetGene[targetGene %in% rownames(numExpr)]
exprSelected <- t(numExpr[matched_genes, , drop = FALSE])

immuneMat <- as.matrix(cibersort_rt)
commonSamples <- intersect(rownames(exprSelected), rownames(immuneMat))
exprSelected <- exprSelected[commonSamples, , drop = FALSE]
immuneMat <- immuneMat[commonSamples, , drop = FALSE]
validImmune <- immuneMat[, apply(immuneMat, 2, sd) > 0, drop = FALSE]

cor_mat <- matrix(NA, nrow = ncol(exprSelected), ncol = ncol(validImmune),
                  dimnames = list(colnames(exprSelected), colnames(validImmune)))
p_val_matrix <- matrix(NA, nrow = ncol(exprSelected), ncol = ncol(validImmune),
                       dimnames = list(colnames(exprSelected), colnames(validImmune)))
for (i in seq_len(ncol(exprSelected))) {
  for (j in seq_len(ncol(validImmune))) {
    test_result <- cor.test(exprSelected[, i], validImmune[, j], method = "spearman", use = "pairwise.complete.obs")
    cor_mat[i, j] <- test_result$estimate
    p_val_matrix[i, j] <- test_result$p.value
  }
}
p_adj_matrix <- matrix(p.adjust(as.vector(p_val_matrix), method = "BH"),
                       nrow = nrow(p_val_matrix), ncol = ncol(p_val_matrix),
                       dimnames = dimnames(p_val_matrix))

corrDF <- data.frame()
for (i in 1:nrow(cor_mat)) {
  for (j in 1:ncol(cor_mat)) {
    corrDF <- rbind(corrDF, data.frame(
      Gene = rownames(cor_mat)[i], ImmuneCell = colnames(cor_mat)[j],
      Correlation = cor_mat[i, j], P.Value = p_val_matrix[i, j],
      P.Adj = p_adj_matrix[i, j],
      Significance = ifelse(p_adj_matrix[i, j] < 0.05, "Yes", "No")))
  }
}
write.csv(corrDF, file = "Gene_Immune_Correlation_Results_CIBERSORT.csv", row.names = FALSE)

bubble_df <- corrDF[!is.na(corrDF$P.Value) & !is.na(corrDF$Correlation), ]
n_genes <- length(unique(bubble_df$Gene))
n_immune <- length(unique(bubble_df$ImmuneCell))
bubble_plot(bubble_df, "ImmuneCell", "ImmuneCell", "Figure6A.pdf",
            max(6, 3 + n_genes * 1.5) - 2, max(8, 2 + n_immune * 0.42) - 4)
immune_abundance <- MCPcounter.estimate(
  expression = numExpr,
  featuresType = "HUGO_symbols",
  genes = read.table(file.path(work_dir, "data", "genes.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
  probesets = read.table(file.path(work_dir, "data", "probesets.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
)
mcp_df <- as.data.frame(t(immune_abundance))
write.csv(mcp_df, file = "Figure6_MCPcounter_Results.csv", row.names = TRUE)

mcpMat <- as.matrix(mcp_df)
commonSamples <- intersect(rownames(exprSelected), rownames(mcpMat))
mcpExpr <- t(numExpr[matched_genes, commonSamples, drop = FALSE])
mcpMat <- mcpMat[commonSamples, , drop = FALSE]
validMcp <- mcpMat[, apply(mcpMat, 2, sd) > 0, drop = FALSE]

mcp_cor <- matrix(NA, nrow = ncol(mcpExpr), ncol = ncol(validMcp),
                  dimnames = list(colnames(mcpExpr), colnames(validMcp)))
mcp_p <- matrix(NA, nrow = ncol(mcpExpr), ncol = ncol(validMcp),
                dimnames = list(colnames(mcpExpr), colnames(validMcp)))
for (i in seq_len(ncol(mcpExpr))) {
  for (j in seq_len(ncol(validMcp))) {
    test_result <- cor.test(mcpExpr[, i], validMcp[, j], method = "spearman", use = "pairwise.complete.obs")
    mcp_cor[i, j] <- test_result$estimate
    mcp_p[i, j] <- test_result$p.value
  }
}
mcp_padj <- matrix(p.adjust(as.vector(mcp_p), method = "BH"),
                   nrow = nrow(mcp_p), ncol = ncol(mcp_p), dimnames = dimnames(mcp_p))

mcpDF <- data.frame()
for (i in 1:nrow(mcp_cor)) {
  for (j in 1:ncol(mcp_cor)) {
    mcpDF <- rbind(mcpDF, data.frame(
      Gene = rownames(mcp_cor)[i], ImmuneCell = colnames(mcp_cor)[j],
      Correlation = mcp_cor[i, j], P.Value = mcp_p[i, j],
      P.Adj = mcp_padj[i, j],
      Significance = ifelse(mcp_padj[i, j] < 0.05, "Yes", "No")))
  }
}
write.csv(mcpDF, file = "Gene_Immune_Correlation_Results_MCPcounter.csv", row.names = FALSE)

mcp_bubble <- mcpDF[!is.na(mcpDF$P.Value) & !is.na(mcpDF$Correlation), ]
n_genes <- length(unique(mcp_bubble$Gene))
n_cell <- length(unique(mcp_bubble$ImmuneCell))
bubble_plot(mcp_bubble, "ImmuneCell", "ImmuneCell", "Figure6B.pdf",
            max(6, 3 + n_genes * 1.5) - 2, max(8, 2 + n_cell * 0.42) - 4)

cytokineGenes <- readLines(cytokine_file)
cytokineGenes <- trimws(cytokineGenes)
cytokineGenes <- cytokineGenes[cytokineGenes != ""]
cytoMatched <- cytokineGenes[cytokineGenes %in% rownames(numExpr)]

hubExpr <- t(numExpr[matched_genes, , drop = FALSE])
cytoExpr <- t(numExpr[cytoMatched, , drop = FALSE])

cyto_cor <- matrix(NA, nrow = length(matched_genes), ncol = length(cytoMatched),
                   dimnames = list(matched_genes, cytoMatched))
cyto_p <- matrix(NA, nrow = length(matched_genes), ncol = length(cytoMatched),
                 dimnames = list(matched_genes, cytoMatched))
for (i in seq_along(matched_genes)) {
  for (j in seq_along(cytoMatched)) {
    test_result <- cor.test(hubExpr[, i], cytoExpr[, j], method = "spearman", use = "pairwise.complete.obs")
    cyto_cor[i, j] <- test_result$estimate
    cyto_p[i, j] <- test_result$p.value
  }
}
cyto_padj <- matrix(p.adjust(as.vector(cyto_p), method = "BH"),
                    nrow = nrow(cyto_p), ncol = ncol(cyto_p), dimnames = dimnames(cyto_p))

cytoDF <- data.frame()
for (i in 1:nrow(cyto_cor)) {
  for (j in 1:ncol(cyto_cor)) {
    cytoDF <- rbind(cytoDF, data.frame(
      Gene = rownames(cyto_cor)[i], Cytokine = colnames(cyto_cor)[j],
      Correlation = cyto_cor[i, j], P.Value = cyto_p[i, j],
      P.Adj = cyto_padj[i, j],
      Significance = ifelse(cyto_padj[i, j] < 0.05, "Yes", "No")))
  }
}
write.csv(cytoDF, file = "Gene_Cytokine_Correlation_Results.csv", row.names = FALSE)

cyto_bubble <- cytoDF[!is.na(cytoDF$P.Value) & !is.na(cytoDF$Correlation), ]
n_genes <- length(unique(cyto_bubble$Gene))
n_cyto <- length(unique(cyto_bubble$Cytokine))
bubble_plot(cyto_bubble, "Cytokine", "Cytokine", "Figure6C.pdf",
            max(6, 3 + n_genes * 1.5) - 2, max(8, 2 + n_cyto * 0.42) - 5)