library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)

work_dir <- getwd()
setwd(work_dir)

gene_exp_file <- "merged_after_batch_removal.csv"
kegg_gmt_file <- file.path(work_dir, "data", "c2.cp.kegg.symbols.gmt")
gene_set_file <- file.path(work_dir, "data", "hub_genes_9.csv")
out_term_num <- 3

expr_df <- read.csv(gene_exp_file, header = TRUE, check.names = FALSE, row.names = 1)
expr_mat <- as.matrix(expr_df)
expr_mat <- matrix(as.numeric(expr_mat), nrow = nrow(expr_mat),
                   dimnames = list(rownames(expr_mat), colnames(expr_mat)))
expr_mat <- limma::avereps(expr_mat)

disease_samples <- grep("_tre$", colnames(expr_mat), value = TRUE)
if (length(disease_samples) == 0) stop("no disease samples ending with _tre")
expr_mat <- expr_mat[, disease_samples, drop = FALSE]

kegg_gmt <- clusterProfiler::read.gmt(kegg_gmt_file)

gene_data <- read.csv(gene_set_file, header = TRUE, stringsAsFactors = FALSE)
gene_list <- as.character(gene_data[[1]])
gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]

plain_num <- function(x) format(x, scientific = FALSE, trim = TRUE)

for (i in seq_along(gene_list)) {
  current_gene <- gene_list[i]
  if (!(current_gene %in% rownames(expr_mat))) next

  this_exp <- expr_mat[current_gene, ]
  median_val <- median(this_exp, na.rm = TRUE)
  if (is.na(median_val)) next
  idx_lo <- which(this_exp < median_val)
  idx_hi <- which(this_exp >= median_val)
  if (length(idx_lo) == 0 || length(idx_hi) == 0) next

  mean_lo <- rowMeans(expr_mat[, idx_lo, drop = FALSE])
  mean_hi <- rowMeans(expr_mat[, idx_hi, drop = FALSE])
  mean_lo[mean_lo < 1e-5] <- 1e-5
  mean_hi[mean_hi < 1e-5] <- 1e-5
  logFC_vec <- log2(mean_hi) - log2(mean_lo)
  logFC_vec <- sort(logFC_vec, decreasing = TRUE)

  gsea_kegg <- tryCatch({
    clusterProfiler::GSEA(logFC_vec, TERM2GENE = kegg_gmt, pvalueCutoff = 1)
  }, error = function(e) NULL)
  if (is.null(gsea_kegg)) next

  kegg_tab <- as.data.frame(gsea_kegg)
  kegg_tab <- kegg_tab[order(-kegg_tab$NES), ]

  out_tab <- kegg_tab
  for (cl in intersect(c("pvalue", "p.adjust", "qvalue"), colnames(out_tab))) {
    out_tab[[cl]] <- plain_num(out_tab[[cl]])
  }
  write.table(out_tab, file.path(work_dir, paste0("Figure5_", current_gene, "_KEGG_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  kegg_up <- head(kegg_tab$ID[kegg_tab$NES > 0], out_term_num)
  kegg_dn <- rev(tail(kegg_tab$ID[kegg_tab$NES < 0], out_term_num))
  kegg_all <- c(kegg_up, kegg_dn)
  if (length(kegg_all) == 0) next

  pal <- brewer.pal(min(length(kegg_all), 8), "Set1")
  if (length(kegg_all) > 8) pal <- rep(pal, length.out = length(kegg_all))

  plt_kegg <- gseaplot2(
    gsea_kegg,
    geneSetID = kegg_all,
    color = pal,
    base_size = 6,
    rel_heights = c(1.2, 0.3, 0.5)
  )
  pdf(file = file.path(work_dir, paste0("Figure5_", current_gene, "_multiterm.pdf")), width = 7.5, height = 5.5)
  print(plt_kegg)
  dev.off()

  if (length(kegg_up) > 0) {
    pdf(file = file.path(work_dir, paste0("Figure5_", current_gene, "_UP.pdf")), width = 7, height = 5.5)
    print(gseaplot2(gsea_kegg, kegg_up, base_size = 12))
    dev.off()
  }
  if (length(kegg_dn) > 0) {
    pdf(file = file.path(work_dir, paste0("Figure5_", current_gene, "_DOWN.pdf")), width = 7, height = 5.5)
    print(gseaplot2(gsea_kegg, kegg_dn, base_size = 12))
    dev.off()
  }
}