library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(limma)
library(ggpubr)
library(rms)
library(regplot)
library(rmda)
library(tools)
library(ggplot2)

work_dir <- getwd()
setwd(work_dir)

source(file.path(work_dir, "data", "ML.R"))

Train_expr <- read.table(file.path(work_dir, "data", "Training_expr.txt"), header = T, sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F)
Train_class <- read.table(file.path(work_dir, "data", "Training_class.txt"), header = T, sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F)
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[, comsam]
Train_class <- Train_class[comsam, , drop = F]

Test_expr <- read.table(file.path(work_dir, "data", "Testing_expr.txt"), header = T, sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F)
Test_class <- read.table(file.path(work_dir, "data", "Testing_class.txt"), header = T, sep = "\t", row.names = 1, check.names = F, stringsAsFactors = F)
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[, comsam]
Test_class <- Test_class[comsam, , drop = F]

comgene <- intersect(rownames(Train_expr), rownames(Test_expr))
Train_expr <- t(Train_expr[comgene, ])
Test_expr <- t(Test_expr[comgene, ])

Train_set <- scaleData(data = Train_expr, centerFlags = T, scaleFlags = T)
Test_set <- scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)

methods <- read.csv(file.path(work_dir, "data", "methods.csv"), stringsAsFactors = FALSE)$Model
methods <- gsub("-| ", "", methods)

classVar <- "outcome"
min.selected.var <- 3

preTrain.var <- list()
set.seed(seed = 777)
Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))
for (method in preTrain.method) {
  preTrain.var[[method]] <- RunML(method = method, Train_set = Train_set,
                                  Train_label = Train_class,
                                  mode = "Variable", classVar = classVar)
}
preTrain.var[["simple"]] <- colnames(Train_set)

model <- list()
set.seed(seed = 777)
Train_set_bk <- Train_set
for (method in methods) {
  method_name <- method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method) == 1) method <- c("simple", method)
  Variable <- preTrain.var[[method[1]]]
  Train_set <- Train_set_bk[, Variable]
  Train_label <- Train_class
  model[[method_name]] <- RunML(method = method[2], Train_set = Train_set,
                                Train_label = Train_label,
                                mode = "Model", classVar = classVar)
  if (length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set <- Train_set_bk
saveRDS(model, "model.rds")

model <- readRDS("model.rds")
methodsValid <- names(model)

RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(fit = model[[method]],
                                       new_data = rbind.data.frame(Train_set, Test_set))
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(data.frame(ID = rownames(RS_mat), RS_mat), "RS_mat.txt", sep = "\t", row.names = F, quote = F)

Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(fit = model[[method]],
                                       new_data = rbind.data.frame(Train_set, Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(data.frame(ID = rownames(Class_mat), Class_mat), "Class_mat.txt", sep = "\t", row.names = F, quote = F)

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit) data.frame(ExtractVar(fit)))
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, "fea_df.txt", sep = "\t", row.names = F, col.names = T, quote = F)

AUC_list <- list()
for (method in methodsValid) {
  AUC_list[[method]] <- RunEval(fit = model[[method]],
                                Test_set = Test_set,
                                Test_label = Test_class,
                                Train_name = "Training",
                                Train_set = Train_set,
                                Train_label = Train_class,
                                cohortVar = "Cohort",
                                classVar = classVar)
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(data.frame(ID = rownames(AUC_mat), AUC_mat), "AUC_mat.txt", sep = "\t", row.names = F, col.names = T, quote = F)

AUC_mat <- read.table("AUC_mat.txt", sep = "\t", row.names = 1, header = T, check.names = F, stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = T)
AUC_mat <- AUC_mat[names(avg_AUC), ]
fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))
if (ncol(AUC_mat) < 3) {
  CohortCol <- c("red", "blue")
} else {
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
}
names(CohortCol) <- colnames(AUC_mat)
cellwidth <- 1
cellheight <- 0.5
hm <- SimpleHeatmap(AUC_mat, avg_AUC, CohortCol, "steelblue",
                    cellwidth = cellwidth, cellheight = cellheight,
                    cluster_columns = F, cluster_rows = F)
pdf("Figure3A.pdf", width = cellwidth * ncol(AUC_mat) + 7, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())
geneFile <- file.path(work_dir, "data", "hub_genes_9.csv")
selectedGenes <- read.csv(geneFile, header = TRUE, stringsAsFactors = FALSE)[[1]]
selectedGenes <- trimws(selectedGenes)
selectedGenes <- selectedGenes[selectedGenes != ""]

test_files <- list.files(pattern = "^test_.*\\.csv$", ignore.case = TRUE)
if (length(test_files) == 0) stop("no test_*.csv files")

my_cols <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666",
  "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
  "#FDB462", "#B3DE69", "#FCCDE5", "#BC80BD", "#CCEBC5"
)

for (test_file in test_files) {
  rt <- read.csv(file.path(work_dir, "data", test_file), header = TRUE, row.names = 1, check.names = FALSE)
  sample_names <- rownames(rt)
  y <- ifelse(grepl("_con", sample_names, ignore.case = TRUE), 0, 1)
  if ("Type" %in% colnames(rt)) rt <- rt[, colnames(rt) != "Type", drop = FALSE]
  if (sum(y == 0) == 0 || sum(y == 1) == 0) next
  found_genes <- selectedGenes[selectedGenes %in% colnames(rt)]
  if (length(found_genes) == 0) next
  rt_filtered <- rt[, found_genes, drop = FALSE]
  while (length(my_cols) < length(found_genes)) {
    my_cols <- c(my_cols, rainbow(length(found_genes) - length(my_cols)))
  }
  dataset_name <- gsub("^test_", "", test_file)
  dataset_name <- gsub("\\.csv$", "", dataset_name, ignore.case = TRUE)
  panel_map <- c(GSE106986 = "B", GSE38974 = "C", GSE8545 = "D")
  panel_letter <- if (dataset_name %in% names(panel_map)) panel_map[[dataset_name]] else "B"
  output_file <- paste0("Figure3", panel_letter, ".pdf")
  pdf(output_file, width = 8, height = 8)
  par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       main = paste("ROC Curves:", dataset_name),
       cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45)
  abline(0, 1, lty = 2, col = "gray70", lwd = 2)
  auc_values <- c()
  for (i in seq_along(found_genes)) {
    gene <- found_genes[i]
    roc_single <- roc(y, as.numeric(rt_filtered[[gene]]), levels = c(0, 1), direction = "auto", quiet = TRUE)
    auc_val <- as.numeric(auc(roc_single))
    auc_values <- c(auc_values, auc_val)
    lines(1 - roc_single$specificities, roc_single$sensitivities, col = my_cols[i], lwd = 2, lty = 1)
  }
  legend_labels <- paste0(found_genes, " (AUC=", sprintf("%.3f", auc_values), ")")
  legend("bottomright", legend = legend_labels, col = my_cols[1:length(found_genes)],
         lwd = 2, lty = 1, bty = "n", cex = 0.85)
  dev.off()
}
input_expr <- file.path(work_dir, "merged_after_batch_removal.csv")
gene_file <- file.path(work_dir, "data", "hub_genes_9.csv")
dat_expr <- read.table(input_expr, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
gene_tbl <- read.csv(gene_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
feature_genes <- as.character(gene_tbl[, 1])
feature_genes <- trimws(feature_genes)
feature_genes <- feature_genes[feature_genes != ""]

dat_expr <- dat_expr[gsub("-", "_", rownames(dat_expr)), , drop = FALSE]
found_genes <- feature_genes %in% rownames(dat_expr)
dat_expr_filt <- dat_expr[feature_genes[found_genes], , drop = FALSE]

df_expr <- as.data.frame(t(dat_expr_filt))
sample_names <- rownames(df_expr)
groups <- gsub("(.*)_([A-Za-z0-9]+)$", "\\2", sample_names)
df_expr$GroupType <- groups

ddinfo <- datadist(df_expr)
options(datadist = "ddinfo")

model_vars <- setdiff(colnames(df_expr), "GroupType")
reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse = " + ")))

lrm_fit <- lrm(reg_formula, data = df_expr, x = TRUE, y = TRUE)

nomo_obj <- nomogram(
  lrm_fit, fun = plogis,
  fun.at = c(0.001, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99), lp = FALSE, funlabel = "Disease Risk"
)
pdf("Figure3E_nomogram.pdf", width = 11, height = 6)
plot(nomo_obj)
dev.off()

median_point <- sapply(df_expr[, model_vars, drop = FALSE], median)
median_df <- as.data.frame(t(median_point))
median_df$GroupType <- 1
rownames(median_df) <- "AllMedian"
regplot(
  lrm_fit, showP = TRUE, rank = "sd", distribution = TRUE,
  observation = median_df, title = "Prediction Nomogram"
)

calibrate_obj <- calibrate(lrm_fit, method = "boot", B = 800)
pdf("Figure3F.pdf", width = 5.5, height = 5.5)
cal_data <- calibrate_obj[, c("predy", "calibrated.orig", "calibrated.corrected")]
plot(cal_data[, "predy"], cal_data[, "calibrated.corrected"],
     type = "l", lwd = 2, col = "black",
     xlab = "Predicted", ylab = "Observed",
     xlim = c(0, 1), ylim = c(0, 1))
abline(a = 0, b = 1, lty = 2, col = "gray50")
lines(cal_data[, "predy"], cal_data[, "calibrated.orig"], lty = 1, lwd = 1, col = "darkblue")
legend("bottomright", legend = c("Ideal", "Apparent", "Bias-corrected"),
       lty = c(2, 1, 1), lwd = c(1, 1, 2), col = c("gray50", "darkblue", "black"), bty = "n")
dev.off()

df_expr$GroupType <- as.numeric(factor(df_expr$GroupType)) - 1
set.seed(123)
df_expr_dca <- df_expr
df_expr_dca$Diagnostic_Model <- predict(lrm_fit, newdata = df_expr, type = "fitted")

dca_combined <- decision_curve(
  GroupType ~ Diagnostic_Model,
  data = df_expr_dca,
  thresholds = seq(0, 1, by = 0.01),
  family = binomial(link = "logit"),
  bootstraps = 100
)

dca_gene_list <- list()
for (g in model_vars) {
  gene_formula <- as.formula(paste("GroupType ~", g))
  dca_gene_list[[g]] <- tryCatch(
    decision_curve(gene_formula, data = df_expr_dca,
                   thresholds = seq(0, 1, by = 0.01),
                   family = binomial(link = "logit"), bootstraps = 100),
    error = function(e) NULL
  )
}
dca_gene_list <- dca_gene_list[!sapply(dca_gene_list, is.null)]
all_dca <- c(list(dca_combined), dca_gene_list)

dca_colors <- c("red", "deepskyblue", "forestgreen", "orange",
                "purple", "magenta", "gold", "brown", "gray40", "cyan")
while (length(dca_colors) < length(all_dca)) {
  dca_colors <- c(dca_colors, rainbow(length(all_dca) - length(dca_colors)))
}
dca_colors <- dca_colors[1:length(all_dca)]
dca_labels <- c("Diagnostic_Model", names(dca_gene_list))

pdf(file = "Figure3G.pdf", width = 9, height = 8)
layout(matrix(c(1, 2), nrow = 2), heights = c(0.8, 5))
par(mar = c(0, 4, 0.5, 2))
plot.new()
legend_labels <- c(dca_labels, "Treat All", "Treat None")
legend_colors <- c(dca_colors, "black", "black")
legend_lty <- c(rep(1, length(dca_labels)), 1, 2)
legend("center", legend = legend_labels, col = legend_colors, lty = legend_lty,
       lwd = 2, ncol = 3, bty = "n", cex = 1.15)
par(mar = c(5, 4, 0.5, 2))
plot_decision_curve(all_dca, xlab = "Threshold Probability", col = dca_colors,
                    confidence.intervals = FALSE, standardize = TRUE,
                    cost.benefit.axis = TRUE, legend.position = "none")
dev.off()