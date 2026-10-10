library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(kernelshap)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(RColorBrewer)
library(pheatmap)
library(gbm)
library(patchwork)

work_dir <- getwd()
setwd(work_dir)

random_seed <- 12345
num_colors <- 20
train_ratio <- 0.7
cv_folds <- 5
pdf_width <- 6
pdf_height <- 6

ml_methods <- data.frame(
  ModelName = c("RF", "XGB", "GBM", "KNN", "LASSO", "Logistic", "NaiveBayes", "C5.0"),
  MethodID = c("rf", "xgbTree", "gbm", "knn", "glmnet", "glm", "naive_bayes", "C5.0"),
  stringsAsFactors = FALSE
)

set.seed(random_seed)

if (num_colors > 9) {
  color_palette <- colorRampPalette(brewer.pal(9, "Pastel1"))(num_colors)
} else {
  color_palette <- brewer.pal(num_colors, "Pastel1")
}

expr_all <- read.csv("merged_after_batch_removal.csv", header = TRUE, check.names = FALSE, row.names = 1)
hub_genes <- read.csv(file.path(work_dir, "data", "hub_genes_9.csv"), header = TRUE, stringsAsFactors = FALSE)[[1]]
hub_genes <- hub_genes[hub_genes %in% rownames(expr_all)]
raw_data <- expr_all[hub_genes, , drop = FALSE]
transposed_data <- t(raw_data)
sample_ids <- rownames(transposed_data)

if (!all(grepl("(_con$|_tre$)", sample_ids, ignore.case = TRUE))) {
  stop("sample names must end with _con or _tre")
}

group_labels <- ifelse(grepl("_con$", sample_ids, ignore.case = TRUE), "Control",
                       ifelse(grepl("_tre$", sample_ids, ignore.case = TRUE), "Treatment", NA))
if (any(is.na(group_labels))) stop("group label failed")

expr_data <- as.data.frame(transposed_data)
expr_data$Group <- as.factor(group_labels)

for (col in colnames(expr_data)[-ncol(expr_data)]) {
  expr_data[[col]] <- as.numeric(expr_data[[col]])
}

expr_data[is.na(expr_data)] <- 0
for (col in colnames(expr_data)[-ncol(expr_data)]) {
  expr_data[[col]][!is.finite(expr_data[[col]])] <- 0
}

feature_vars <- apply(expr_data[, -ncol(expr_data)], 2, var, na.rm = TRUE)
zero_var_features <- names(feature_vars[feature_vars == 0 | is.na(feature_vars)])
if (length(zero_var_features) > 0) {
  expr_data <- expr_data[, !(colnames(expr_data) %in% zero_var_features)]
}

if (length(unique(expr_data$Group)) != 2) stop("group number is not 2")
split_idx <- createDataPartition(y = expr_data$Group, p = train_ratio, list = FALSE)
train_set <- expr_data[split_idx, ]
test_set <- expr_data[-split_idx, ]
test_labels <- test_set$Group
test_features <- test_set[, -ncol(test_set)]

auc_results <- c()
model_auc_map <- list()
roc_colors <- color_palette[1:nrow(ml_methods)]
roc_ci_list <- list()
trained_models <- list()
first_plot <- TRUE

pdf(file = "Figure4A.pdf", width = pdf_width, height = pdf_height)
par(mar = c(5, 5, 4, 2) + 0.1)
for (i in 1:nrow(ml_methods)) {
  mdl_name <- ml_methods$ModelName[i]
  mdl_id <- ml_methods$MethodID[i]
  if (nrow(train_set) < 10) stop("training set too small")

  mdl <- tryCatch({
    if (mdl_id == "xgbTree") {
      caret::train(Group ~ ., data = train_set, method = mdl_id,
                   trControl = caret::trainControl(method = "cv", number = cv_folds,
                                                   savePredictions = TRUE, allowParallel = FALSE,
                                                   verboseIter = FALSE),
                   tuneGrid = expand.grid(nrounds = c(50, 100),
                                          max_depth = c(3, 6),
                                          eta = 0.3,
                                          gamma = 0,
                                          colsample_bytree = 1,
                                          min_child_weight = 1,
                                          subsample = 1),
                   verbose = FALSE)
    } else {
      caret::train(Group ~ ., data = train_set, method = mdl_id,
                   trControl = caret::trainControl(method = "repeatedcv", number = cv_folds,
                                                   savePredictions = TRUE, allowParallel = FALSE))
    }
  }, error = function(e) NULL)

  if (is.null(mdl)) next
  trained_models[[mdl_id]] <- mdl

  pred_prob <- predict(mdl, newdata = test_features, type = "prob")
  if (!"Treatment" %in% colnames(pred_prob)) next
  roc_obj <- roc(as.numeric(test_labels) - 1, as.numeric(pred_prob[, "Treatment"]), quiet = TRUE)
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj, conf.level = 0.95)
  roc_ci_list[[mdl_id]] <- auc_ci
  auc_results <- c(auc_results, sprintf("%s: %.03f [%.03f, %.03f]", mdl_name, auc_val, auc_ci[1], auc_ci[3]))
  model_auc_map[[mdl_id]] <- auc_val
  if (first_plot) {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "", col = roc_colors[i], lwd = 3)
    first_plot <- FALSE
  } else {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "", col = roc_colors[i], lwd = 3, add = TRUE)
  }
}
legend("bottomright", auc_results, col = roc_colors[1:length(auc_results)], lwd = 3, bty = "n", cex = 0.9)
dev.off()
auc_vec <- unlist(model_auc_map)
if (length(auc_vec) == 0) stop("AUC result is empty")
best_method <- names(which.max(auc_vec))
final_model <- caret::train(Group ~ ., data = train_set, method = best_method,
                            trControl = caret::trainControl(method = "repeatedcv", number = cv_folds,
                                                            savePredictions = TRUE))
shap_fit <- kernelshap(
  object = final_model,
  X = train_set[, -ncol(train_set)],
  pred_fun = function(obj, newdata) {
    pred <- predict(obj, newdata, type = "prob")
    if (!"Treatment" %in% colnames(pred)) stop("no Treatment column in prediction")
    pred[, "Treatment"]
  }
)

shap_obj <- shapviz(shap_fit, X_pred = train_set[, -ncol(train_set)],
                    X = train_set[, -ncol(train_set)], interactions = TRUE)
shap_importance <- sort(colMeans(abs(shap_obj$S)), decreasing = TRUE)
top_features <- names(shap_importance)

custom_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "right")

options(shapviz.colors = color_palette)

pdf(file = "Figure4B.pdf", width = pdf_width, height = pdf_height)
print(sv_importance(shap_obj, kind = "bar", show_numbers = TRUE) + custom_theme)
dev.off()

pdf(file = "Figure4C.pdf", width = pdf_width, height = pdf_height)
print(sv_importance(shap_obj, kind = "bee", show_numbers = TRUE) + custom_theme)
dev.off()

dependence_plots_list <- list()
for (i in seq_along(top_features)) {
  gene <- top_features[i]
  dependence_plot <- sv_dependence(shap_obj, v = gene) +
    custom_theme +
    ggtitle(gene) +
    theme(plot.title = element_text(size = 11, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_blank())
  dependence_plots_list[[i]] <- dependence_plot
}

n_cols <- 4
n_rows <- ceiling(length(dependence_plots_list) / n_cols)
pdf(file = "Figure4D.pdf", width = 16, height = 4 * n_rows)
print(wrap_plots(dependence_plots_list, ncol = n_cols))
dev.off()

pdf(file = "Figure4E.pdf", width = pdf_width, height = pdf_height)
print(sv_waterfall(shap_obj, row_id = 30) + custom_theme)
dev.off()

pdf(file = "Figure4F.pdf", width = pdf_width, height = pdf_height)
print(sv_force(shap_obj, row_id = 30) + custom_theme)
dev.off()

top20_genes <- head(shap_importance, 20)
write.csv(data.frame(Gene = names(top20_genes), Importance = as.numeric(top20_genes)),
          "shap_top_genes.csv", row.names = FALSE)