workDir <- getwd()

library(Seurat)
library(harmony)
library(limma)
library(SingleR)
library(celldex)

library(DoubletFinder)
library(celda)
library(dplyr)
library(magrittr)
library(tidyr)
library(stringr)

library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(viridis)
library(scales)
library(patchwork)

library(rstatix)

library(DT)
library(progress)

library(openxlsx)

library(ComplexHeatmap)
library(ggsci)

library(Matrix)
library(svglite)
library(clustree)
library(cluster)
library(ggrepel)
library(ggunchull)
library(scRNAtoolVis)

library(plot1cell)
library(circlize)

message("All required R packages loaded successfully!")

raw_10x_dir <- file.path(workDir, "data", "raw_10x")
sample_dirs <- list.dirs(raw_10x_dir, recursive = FALSE)
if (length(sample_dirs) == 0) {
  stop("No 10X sample directories found in: ", raw_10x_dir)
}

counts_list <- list()
for (i in seq_along(sample_dirs)) {
  sample_name <- basename(sample_dirs[i])
  if (grepl("^[Cc]ontrol", sample_name)) sample_name <- "Control"
  if (grepl("^[Dd]isease", sample_name)) sample_name <- "Disease"
  mat <- tryCatch(Read10X(data.dir = sample_dirs[i]), error = function(e) NULL)
  if (is.null(mat)) next
  colnames(mat) <- paste0(sample_name, ".", colnames(mat))
  counts_list[[sample_name]] <- mat
}
if (length(counts_list) == 0) stop("No 10X data could be read")
counts_combined <- do.call(cbind, counts_list)

sc_obj <- CreateSeuratObject(counts_combined, min.cells = 10, min.features = 40)
single_cell_data <- GetAssayData(sc_obj, slot = "counts")
saveRDS(single_cell_data, file = file.path(workDir, "scRNA_matrix.rds"))
message(sprintf("scRNA_matrix.rds written: %d genes x %d cells",
                nrow(single_cell_data), ncol(single_cell_data)))

input_file <- file.path(workDir, "scRNA_matrix.rds")
showGenes <- read.csv(file.path(workDir, "data", "hub_genes_9.csv"), stringsAsFactors = FALSE)[, 1]
showGenes <- trimws(showGenes)
showGenes <- showGenes[showGenes != ""]

analysis_params <- list(
  logFCfilter = 1,
  adjPvalFilter = 0.05,
  min_cells = 10,
  min_features = 30,

  nFeature_RNA_min = 20,
  percent_mt_max = 20,

  normalization_method = "LogNormalize",
  scale_factor = 10000,
  selection_method = "vst",
  n_variable_features = 1500,

  n_pcs = 20,
  pcSelect = 13,
  cluster_resolution = 0.2,
  cluster_resolution_range = seq(0.1, 1.0, 0.1),

  jackstraw_replicates = 100,
  cumulative_variance_threshold = 0.90,
  small_cluster_threshold = 50,
  max_clusters_silhouette = 20
)

message("Reading single cell data...")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

message(sprintf("Reading data from: %s", basename(input_file)))
single_cell_data <- readRDS(input_file)

message(sprintf("Data loaded successfully: %d genes x %d cells",
                nrow(single_cell_data), ncol(single_cell_data)))

control_cells_found <- any(grepl("^Control\\.", colnames(single_cell_data), ignore.case = FALSE))
disease_cells_found <- any(grepl("^Disease\\.", colnames(single_cell_data), ignore.case = FALSE))

if(control_cells_found) {
  message("Detected Control prefix in cells")
}
if(disease_cells_found) {
  message("Detected Disease prefix in cells")
}

if(!control_cells_found && !disease_cells_found) {
  stop("ERROR: No Control or Disease prefixes found in cell names! Please ensure your data has cells with 'Control.' or 'Disease.' prefixes")
}

if(!control_cells_found) {
  stop("ERROR: No Control group cells found!")
}
if(!disease_cells_found) {
  stop("ERROR: No Disease group cells found!")
}

if (ncol(single_cell_data) == 0) {
  stop("Data file is empty!")
}

message(sprintf("Data contains %d genes and %d cells", nrow(single_cell_data), ncol(single_cell_data)))


control_cells <- sum(grepl("^Control\\.", colnames(single_cell_data)))
disease_cells <- sum(grepl("^Disease\\.", colnames(single_cell_data)))
message(sprintf("Control cells: %d, Disease cells: %d", control_cells, disease_cells))


saveRDS(single_cell_data, file.path(workDir, "Raw_Combined_Data.rds"))

message("Creating Seurat object and performing quality control.")
message("Seurat...")

Peripheral_Blood_Mononuclear_Cells <- CreateSeuratObject(
  counts = single_cell_data,
  min.cells = analysis_params$min_cells,
  min.features = analysis_params$min_features
)

message("Seurat object created")

Peripheral_Blood_Mononuclear_Cells$Type <- ifelse(
  grepl("^Control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)), "Control", "Disease"
)
Peripheral_Blood_Mononuclear_Cells[["percent.mt"]] <- PercentageFeatureSet(
  object = Peripheral_Blood_Mononuclear_Cells,
  pattern = "^MT-"
)

pdf(file = file.path(workDir, "FigureS1A_QC_Before_Filter.pdf"),
    width = 18, height = 6.5)
print(
  VlnPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
  )
)
dev.off()

cells_before <- ncol(Peripheral_Blood_Mononuclear_Cells)
genes_before <- nrow(Peripheral_Blood_Mononuclear_Cells)

Peripheral_Blood_Mononuclear_Cells <- subset(
  x = Peripheral_Blood_Mononuclear_Cells,
  subset = nFeature_RNA > analysis_params$nFeature_RNA_min &
           percent.mt < analysis_params$percent_mt_max
)

pdf(file = file.path(workDir, "FigureS1A_QC_After_Filter.pdf"),
    width = 18, height = 6.5)
p_after <- VlnPlot(
  object = Peripheral_Blood_Mononuclear_Cells,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)
p_after[[3]] <- p_after[[3]] + ylim(0, 100)
print(p_after)
dev.off()

cells_after <- ncol(Peripheral_Blood_Mononuclear_Cells)
genes_after <- nrow(Peripheral_Blood_Mononuclear_Cells)

message(sprintf("Before filtering: %d cells, %d genes", cells_before, genes_before))
message(sprintf("Before QC: %d cells / %d genes", cells_before, genes_before))
message(sprintf("After filtering: %d cells, %d genes", cells_after, genes_after))
message(sprintf("After QC: %d cells / %d genes", cells_after, genes_after))
message(sprintf("Retention rate: cells %.2f%%, genes %.2f%%",
                (cells_after/cells_before)*100, (genes_after/genes_before)*100))
message(sprintf("Retained: %.2f%% cells / %.2f%% genes",
                (cells_after/cells_before)*100, (genes_after/genes_before)*100))

qc_summary <- data.frame(
  Metric = c("Cells_Before_Filter", "Cells_After_Filter", "Genes_Before_Filter",
             "Genes_After_Filter", "Cell_Retention_Rate", "Gene_Retention_Rate"),
  Value = c(cells_before, cells_after, genes_before, genes_after,
            paste0(round((cells_after/cells_before)*100, 2), "%"),
            paste0(round((genes_after/genes_before)*100, 2), "%"))
)

write.csv(qc_summary, file.path(workDir, "QC_Summary.csv"),
          row.names = FALSE)

message("decontXRNA...")

raw_counts <- GetAssayData(Peripheral_Blood_Mononuclear_Cells, layer = "counts")

decontX_results <- decontX(raw_counts)

decontXcounts <- decontX_results$decontXcounts
contamination <- decontX_results$contamination

Peripheral_Blood_Mononuclear_Cells$decontX_contamination <- contamination

message(sprintf("Mean ambient RNA contamination: %.2f%%", mean(contamination) * 100))
message(sprintf("Median ambient RNA contamination: %.2f%%", median(contamination) * 100))
message(sprintf("Maximum ambient RNA contamination: %.2f%%", max(contamination) * 100))

pdf(file.path(workDir, "FigureS1A_decontX_Contamination.pdf"), width = 8, height = 6)
hist(contamination, breaks = 50, main = "decontX Contamination Distribution",
     xlab = "Contamination Proportion", col = "steelblue", border = "white")
abline(v = mean(contamination), col = "red", lwd = 2, lty = 2)
legend("topright", legend = paste("Mean:", round(mean(contamination), 3)), col = "red", lty = 2, lwd = 2)
dev.off()

Peripheral_Blood_Mononuclear_Cells[["RNA"]]$counts <- decontXcounts

decontX_summary <- data.frame(
  Metric = c("Mean_Contamination", "Median_Contamination", "Max_Contamination", "Min_Contamination"),
  Value = c(round(mean(contamination), 4), round(median(contamination), 4),
            round(max(contamination), 4), round(min(contamination), 4))
)
write.csv(decontX_summary, file.path(workDir, "decontX_Summary.csv"), row.names = FALSE)

high_contamination_threshold <- 0.5
cells_before_contam_filter <- ncol(Peripheral_Blood_Mononuclear_Cells)

Peripheral_Blood_Mononuclear_Cells <- subset(
  Peripheral_Blood_Mononuclear_Cells,
  subset = decontX_contamination < high_contamination_threshold
)

cells_after_contam_filter <- ncol(Peripheral_Blood_Mononuclear_Cells)
message(sprintf("Cells with contamination > 50%%: %d",
                cells_before_contam_filter - cells_after_contam_filter))
message(sprintf("Cells retained after decontX: %d", cells_after_contam_filter))

message("decontXRNA")

message("DoubletFinder...")

cells_before_doublet <- ncol(Peripheral_Blood_Mononuclear_Cells)

tmp_group <- ifelse(grepl("^Control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)), "Control", "Disease")
Peripheral_Blood_Mononuclear_Cells$tmp_group <- tmp_group

run_doubletfinder_batch <- function(seurat_obj) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

  sweep.res <- paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_val <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  n <- ncol(seurat_obj)
  dr <- n * 0.8 / 100000
  hp <- modelHomotypic(seurat_obj$seurat_clusters)
  nExp <- round(round(dr * n) * (1 - hp))

  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = pK_val, nExp = nExp, sct = FALSE)

  df_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)[1]
  seurat_obj$Doublet_Classification <- seurat_obj@meta.data[[df_col]]
  return(seurat_obj)
}

doublet_results <- list()
for (grp in c("Control", "Disease")) {
  message(sprintf("   %s ...", grp))
  grp_obj <- subset(Peripheral_Blood_Mononuclear_Cells, subset = tmp_group == grp)
  message(sprintf("  %s : %d", grp, ncol(grp_obj)))

  tryCatch({
    grp_obj <- run_doubletfinder_batch(grp_obj)
    doublet_results[[grp]] <- grp_obj$Doublet_Classification
    names(doublet_results[[grp]]) <- colnames(grp_obj)
    n_doublet <- sum(grp_obj$Doublet_Classification == "Doublet")
    message(sprintf("  %s  %d ", grp, n_doublet))
  }, error = function(e) {
    stop(sprintf("%s DoubletFinder: %s", grp, e$message))
  })
}

all_doublet_class <- c(doublet_results[["Control"]], doublet_results[["Disease"]])
Peripheral_Blood_Mononuclear_Cells$Doublet_Classification <- all_doublet_class[colnames(Peripheral_Blood_Mononuclear_Cells)]

Peripheral_Blood_Mononuclear_Cells <- NormalizeData(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(Peripheral_Blood_Mononuclear_Cells, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- RunPCA(Peripheral_Blood_Mononuclear_Cells, verbose = FALSE)
Peripheral_Blood_Mononuclear_Cells <- RunUMAP(Peripheral_Blood_Mononuclear_Cells, dims = 1:10, verbose = FALSE)

pdf(file.path(workDir, "FigureS1A_DoubletFinder_UMAP.pdf"), width = 10, height = 8)
print(DimPlot(Peripheral_Blood_Mononuclear_Cells, group.by = "Doublet_Classification", cols = c("Singlet" = "grey", "Doublet" = "red")))
dev.off()

Peripheral_Blood_Mononuclear_Cells <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Doublet_Classification == "Singlet")

Peripheral_Blood_Mononuclear_Cells$tmp_group <- NULL

cells_after_doublet <- ncol(Peripheral_Blood_Mononuclear_Cells)
message(sprintf("Cells before doublet removal: %d", cells_before_doublet))
message(sprintf("Cells after doublet removal: %d", cells_after_doublet))
message(sprintf("Doublets removed: %d (%.2f%%)", cells_before_doublet - cells_after_doublet, (cells_before_doublet - cells_after_doublet)/cells_before_doublet*100))

Peripheral_Blood_Mononuclear_Cells@reductions <- list()
Peripheral_Blood_Mononuclear_Cells@graphs <- list()
Peripheral_Blood_Mononuclear_Cells$seurat_clusters <- NULL

message("Starting data preprocessing")

Peripheral_Blood_Mononuclear_Cells <- NormalizeData(
  object = Peripheral_Blood_Mononuclear_Cells,
  normalization.method = analysis_params$normalization_method,
  scale.factor = analysis_params$scale_factor
)
message("Data normalization completed")

Peripheral_Blood_Mononuclear_Cells <- FindVariableFeatures(
  object = Peripheral_Blood_Mononuclear_Cells,
  selection.method = analysis_params$selection_method,
  nfeatures = analysis_params$n_variable_features
)
message(sprintf("Identified %d variable features", length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells))))
message(sprintf(" %d ", length(VariableFeatures(Peripheral_Blood_Mononuclear_Cells))))

top10_genes <- head(VariableFeatures(Peripheral_Blood_Mononuclear_Cells), 10)
message("Top 10 variable genes")

pdf(file = file.path(workDir, "FigureS1A_Variable_Features.pdf"),
    width = 15, height = 8)
plot1 <- VariableFeaturePlot(Peripheral_Blood_Mononuclear_Cells)
plot2 <- LabelPoints(plot = plot1, points = top10_genes, repel = TRUE)
print(plot1 + plot2)
dev.off()

write.csv(
  data.frame(Gene = VariableFeatures(Peripheral_Blood_Mononuclear_Cells)),
  file.path(workDir, "Variable_Genes_List.csv"),
  row.names = FALSE
)

Peripheral_Blood_Mononuclear_Cells <- ScaleData(Peripheral_Blood_Mononuclear_Cells)
message("Data scaling completed")

message("Performing PCA analysis")

Peripheral_Blood_Mononuclear_Cells <- RunPCA(
  object = Peripheral_Blood_Mononuclear_Cells,
  npcs = analysis_params$n_pcs,
  features = VariableFeatures(object = Peripheral_Blood_Mononuclear_Cells)
)

message("Performing Harmony batch correction using orig.ident (8 samples) as batch...")

Peripheral_Blood_Mononuclear_Cells <- RunHarmony(
  object = Peripheral_Blood_Mononuclear_Cells,
  group.by.vars = "orig.ident",
  reduction.use = "pca",
  dims.use = 1:analysis_params$n_pcs,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)

write.csv(
  Embeddings(Peripheral_Blood_Mononuclear_Cells, "harmony"),
  file = file.path(workDir, "Harmony_Embedding.csv")
)
message("Harmony embedding saved to 05_Harmony_Embedding.csv")

p_pca_before <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "pca",
  group.by = "orig.ident",
  raster = FALSE
) + ggtitle("PCA before batch correction (by sample)")

p_harmony_after <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "harmony",
  group.by = "orig.ident",
  raster = FALSE
) + ggtitle("Harmony after batch correction (by sample)")

pdf(file = file.path(workDir, "FigureS1B.pdf"),
    width = 14, height = 6)
print(p_pca_before | p_harmony_after)
dev.off()
message("PCA vs Harmony comparison (by sample) saved.")

plot_group_tmp <- ifelse(
  grepl("^Control\\.", colnames(Peripheral_Blood_Mononuclear_Cells)),
  "Control", "Disease"
)
Peripheral_Blood_Mononuclear_Cells$Group_plot_tmp <- plot_group_tmp

p_pca_before_grp <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "pca",
  group.by = "Group_plot_tmp",
  raster = FALSE
) + ggtitle("PCA before correction (by group)")

p_harmony_after_grp <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "harmony",
  group.by = "Group_plot_tmp",
  raster = FALSE
) + ggtitle("Harmony after correction (by group)")

pdf(file = file.path(workDir, "FigureS1B_byGroup.pdf"),
    width = 14, height = 6)
print(p_pca_before_grp | p_harmony_after_grp)
dev.off()
message("PCA vs Harmony comparison (by group) saved.")

Peripheral_Blood_Mononuclear_Cells$Group_plot_tmp <- NULL

pdf(file = file.path(workDir, "FigureS1B_PCA_Plot.pdf"),
    width = 8, height = 6)
print(DimPlot(Peripheral_Blood_Mononuclear_Cells, reduction = "pca") +
      NoGrid())
dev.off()

for (i in 1:4) {
  pdf(file = file.path(workDir, paste0("FigureS1C_PCA_Heatmap_PC", i, ".pdf")),
      width = 10, height = 8)
  print(
    DimHeatmap(
      Peripheral_Blood_Mononuclear_Cells,
      dims = i,
      cells = 500,
      balanced = TRUE
    )
  )
  dev.off()
}

Peripheral_Blood_Mononuclear_Cells <- JackStraw(
  Peripheral_Blood_Mononuclear_Cells,
  num.replicate = analysis_params$jackstraw_replicates
)
Peripheral_Blood_Mononuclear_Cells <- ScoreJackStraw(
  Peripheral_Blood_Mononuclear_Cells,
  dims = 1:analysis_params$n_pcs
)

pdf(file = file.path(workDir, "FigureS1C.pdf"),
    width = 9, height = 6)
print(JackStrawPlot(Peripheral_Blood_Mononuclear_Cells, dims = 1:15))
dev.off()

pca_results <- Peripheral_Blood_Mononuclear_Cells[["pca"]]@stdev
pca_variance <- pca_results^2
cumulative_variance <- cumsum(pca_variance)

pca_table <- data.frame(
  PC = 1:length(pca_variance),
  StdDev = pca_results,
  Variance = pca_variance,
  CumulativeVariance = cumulative_variance
)

write.csv(pca_table, file.path(workDir, "PCA_Results.csv"),
          row.names = FALSE)

selected_pcs <- which(cumulative_variance >= analysis_params$cumulative_variance_threshold)[1]
if(is.na(selected_pcs)) selected_pcs <- analysis_params$pcSelect
message(sprintf("Selected %d principal components based on cumulative variance", selected_pcs))
message(sprintf(" %d ", selected_pcs))

message("Performing clustering analysis")

Peripheral_Blood_Mononuclear_Cells <- FindNeighbors(
  object = Peripheral_Blood_Mononuclear_Cells,
  reduction = "harmony",
  dims = 1:analysis_params$pcSelect
)

message("Performing multi-resolution clustering")
for (res in analysis_params$cluster_resolution_range) {
  Peripheral_Blood_Mononuclear_Cells <- FindClusters(
    Peripheral_Blood_Mononuclear_Cells,
    resolution = res,
    verbose = FALSE
  )
}

pdf(file.path(workDir, "FigureS1C_Clustree.pdf"),
    width = 15, height = 10)
clustree_plot <- clustree(Peripheral_Blood_Mononuclear_Cells@meta.data, prefix = "RNA_snn_res.") +
  theme_classic() +
  labs(title = "Clustering Tree Across Different Resolutions")
print(clustree_plot)
dev.off()

Peripheral_Blood_Mononuclear_Cells <- FindClusters(
  object = Peripheral_Blood_Mononuclear_Cells,
  resolution = analysis_params$cluster_resolution
)
message(sprintf("Clustering completed with resolution %.1f, identified %d clusters",
                analysis_params$cluster_resolution,
                length(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))))
message(sprintf(" %.1f  %d ",
                analysis_params$cluster_resolution,
                length(unique(Peripheral_Blood_Mononuclear_Cells$seurat_clusters))))

Peripheral_Blood_Mononuclear_Cells <- RunUMAP(
  object = Peripheral_Blood_Mononuclear_Cells,
  reduction = "pca",
  dims = 1:analysis_params$pcSelect,
  reduction.name = "umap.pca"
)

Peripheral_Blood_Mononuclear_Cells <- RunUMAP(
  object = Peripheral_Blood_Mononuclear_Cells,
  reduction = "harmony",
  dims = 1:analysis_params$pcSelect,
  reduction.name = "umap"
)

p_umap_before <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "umap.pca",
  group.by = "orig.ident",
  raster = FALSE
) + ggtitle("UMAP before batch correction (by sample)")

p_umap_after <- DimPlot(
  Peripheral_Blood_Mononuclear_Cells,
  reduction = "umap",
  group.by = "orig.ident",
  raster = FALSE
) + ggtitle("UMAP after Harmony correction (by sample)")

pdf(file = file.path(workDir, "FigureS1B_UMAP_Before_vs_After_Harmony.pdf"),
    width = 14, height = 6)
print(p_umap_before | p_umap_after)
dev.off()
message("UMAP before/after Harmony comparison saved.")

message("Generating UMAP visualization")

pdf(file = file.path(workDir, "Figure9A.pdf"),
    width = 10, height = 8)
print(
  DimPlot(
    object = Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    label.size = 4
  ) +
  ggtitle("UMAP Clustering Results") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()

write.table(
  Peripheral_Blood_Mononuclear_Cells$seurat_clusters,
  file = file.path(workDir, "Cluster_Assignments.txt"),
  quote = FALSE,
  sep = "\t",
  col.names = FALSE
)

message("Performing automatic cell type annotation using SingleR...")
message("SingleR...")

if (!requireNamespace("SingleR", quietly = TRUE) || !requireNamespace("celldex", quietly = TRUE)) {
  message("Warning: SingleR or celldex package not installed, using default annotation")
  message("SingleR  celldex ")

  default_labels <- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
  Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- default_labels
  Peripheral_Blood_Mononuclear_Cells$cell_type <- default_labels
  Idents(Peripheral_Blood_Mononuclear_Cells) <- default_labels

} else {

  tryCatch({
    message("Extracting expression data and cluster information...")

    counts <- GetAssayData(object = Peripheral_Blood_Mononuclear_Cells, layer = "data")
    clusters <- Peripheral_Blood_Mononuclear_Cells@meta.data$seurat_clusters

    message("Expression matrix dimensions: ", nrow(counts), " genes x ", ncol(counts), " cells")
    message("Reference matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")
    message("Number of clusters: ", length(unique(clusters)))
    message("Number of clusters: ", length(unique(clusters)))

    message("Loading HumanPrimaryCellAtlasData reference...")
    message("HumanPrimaryCellAtlasData...")
    ref <- celldex::HumanPrimaryCellAtlasData()

    message("Reference dataset loaded: ", nrow(ref), " genes x ", ncol(ref), " samples")
    message("Reference: ", nrow(ref), " x ", ncol(ref))

    message("Running SingleR cluster-level annotation...")
    message("SingleR...")

    singler_clusters <- SingleR(
      test = counts,
      ref = ref,
      labels = ref$label.main,
      clusters = clusters
    )

    clusterAnn <- data.frame(
      Cluster = rownames(singler_clusters),
      CellType = singler_clusters$labels,
      stringsAsFactors = FALSE
    )

    write.csv(clusterAnn,
              file = file.path(workDir, "Cluster_to_CellType_Mapping.csv"),
              row.names = FALSE)

    message("Cluster annotation completed. Results:")
    for(i in 1:nrow(clusterAnn)) {
      message(sprintf("  Cluster %s -> %s", clusterAnn$Cluster[i], clusterAnn$CellType[i]))
    }

    message("Mapping cluster annotations to individual cells...")

    cluster_to_celltype <- setNames(clusterAnn$CellType, clusterAnn$Cluster)

    cell_clusters <- Peripheral_Blood_Mononuclear_Cells@meta.data$seurat_clusters

    cell_types_from_clusters <- cluster_to_celltype[as.character(cell_clusters)]
    names(cell_types_from_clusters) <- colnames(Peripheral_Blood_Mononuclear_Cells)

    if(any(is.na(cell_types_from_clusters))) {
      message("Warning: Some cells have missing annotations, using cluster names as fallback")
      cell_types_from_clusters[is.na(cell_types_from_clusters)] <- paste0("Cluster_",
        cell_clusters[is.na(cell_types_from_clusters)])
    }

    message("Adding cell type annotations to Seurat object...")
    message("Seurat...")

    Peripheral_Blood_Mononuclear_Cells <- AddMetaData(
      Peripheral_Blood_Mononuclear_Cells,
      metadata = cell_types_from_clusters,
      col.name = "cell_type"
    )

    Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- cell_types_from_clusters

    Idents(Peripheral_Blood_Mononuclear_Cells) <- cell_types_from_clusters

    message("Verification:")
    message("  SingleR_labels column exists: ", "SingleR_labels" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  cell_type column exists: ", "cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  Number of unique cell types: ", length(unique(cell_types_from_clusters)))
    message("  : ", length(unique(cell_types_from_clusters)))

    cell_type_table <- table(cell_types_from_clusters)
    message("\nCell type distribution:")
    for(ct_name in names(cell_type_table)) {
      message(sprintf("  %s: %d cells", ct_name, cell_type_table[ct_name]))
    }

    message("SingleR annotation completed successfully!")
    message("SingleR")

  }, error = function(e) {
    message("ERROR: SingleR annotation failed, using default annotation")
    message(sprintf("SingleR: %s", e$message))
    message("Error details: ", conditionMessage(e))
    message("Error: ", conditionMessage(e))

    message("Creating default cluster-based annotations...")

    default_labels <- paste0("Cluster_", Peripheral_Blood_Mononuclear_Cells$seurat_clusters)
    Peripheral_Blood_Mononuclear_Cells$SingleR_labels <- default_labels
    Peripheral_Blood_Mononuclear_Cells$cell_type <- default_labels
    Idents(Peripheral_Blood_Mononuclear_Cells) <- default_labels

    message("Default annotation verification:")
    message("  SingleR_labels column exists: ", "SingleR_labels" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
    message("  cell_type column exists: ", "cell_type" %in% colnames(Peripheral_Blood_Mononuclear_Cells@meta.data))
  })
}

message(sprintf("Automatic annotation completed, identified %d cell types",
                length(unique(Idents(Peripheral_Blood_Mononuclear_Cells)))))
message(sprintf(" %d ",
                length(unique(Idents(Peripheral_Blood_Mononuclear_Cells)))))

cell_type_counts <- table(Idents(Peripheral_Blood_Mononuclear_Cells))
message("\nCell type statistics")
for(ct in names(sort(cell_type_counts, decreasing = TRUE))) {
  message(sprintf("  %s: %d cells", ct, cell_type_counts[ct]))
}

pdf(file = file.path(workDir, "Figure9B.pdf"),
    width = 12, height = 8)
print(
  DimPlot(
    Peripheral_Blood_Mononuclear_Cells,
    reduction = "umap",
    pt.size = 1.5,
    label = TRUE,
    repel = TRUE
  ) +
  ggtitle("Automatic Cell Type Annotation (SingleR)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold")
  )
)
dev.off()

message("Performing cell proportion analysis")

cell_type_counts <- table(Idents(Peripheral_Blood_Mononuclear_Cells),
                          Peripheral_Blood_Mononuclear_Cells$Type)

prop_within_group <- prop.table(cell_type_counts, margin = 2) * 100

proportion_data <- as.data.frame(cell_type_counts)
colnames(proportion_data) <- c("CellType", "Group", "Count")

proportion_summary <- proportion_data %>%
  dplyr::group_by(Group) %>%
  dplyr::mutate(
    Total = sum(Count),
    Proportion = Count / Total
  ) %>%
  dplyr::ungroup()

write.csv(proportion_summary,
          file.path(workDir, "Cell_Type_Proportions.csv"),
          row.names = FALSE)

pdf(file.path(workDir, "Figure9C.pdf"),
    width = 8, height = 6)

unique_celltypes <- unique(proportion_summary$CellType)
nColors <- length(unique_celltypes)
celltype_colors <- colorRampPalette(brewer.pal(min(12, nColors), "Set3"))(nColors)
names(celltype_colors) <- unique_celltypes

p <- ggplot(proportion_summary, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.3) +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  labs(title = "Cell Type Proportions: Control vs Disease", y = "Proportion")

print(p)
dev.off()

message("Generating expression plots for genes of interest...")

available_genes <- showGenes[showGenes %in% rownames(GetAssayData(Peripheral_Blood_Mononuclear_Cells, assay = "RNA", layer = "data"))]
missing_genes <- setdiff(showGenes, available_genes)

if (length(missing_genes) > 0) {
  message(paste("The following genes were not found in the dataset:", paste(missing_genes, collapse = ", ")))
  message(paste(":", paste(missing_genes, collapse = ", ")))
}

if (length(available_genes) > 0) {
  pdf(file = file.path(workDir, "Supp_ShowGenes_Violin.pdf"),
      width = 15, height = 10)
  print(
    VlnPlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes,
      ncol = min(3, length(available_genes))
    ) +
    ggtitle("Gene Expression Distribution Across Cell Types")
  )
  dev.off()

  pdf(file = file.path(workDir, "Supp_ShowGenes_Feature.pdf"),
      width = 15, height = 12)
  print(
    FeaturePlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes,
      cols = c("#E5E5E5", "#FF0000"),
      reduction = "umap",
      ncol = min(3, length(available_genes))
    ) & theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold")
    )
  )
  dev.off()

  pdf(file = file.path(workDir, "Figure9E.pdf"),
      width = 10, height = 8)
  print(
    DotPlot(
      object = Peripheral_Blood_Mononuclear_Cells,
      features = available_genes
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Expression of Interest Genes Across Cell Types",
         x = "Genes", y = "Cell Types")
  )
  dev.off()
}

message("\n================================================================================")
message("         Step 18: Single Cell Gene Knockout Simulation Analysis                 ")


if (!requireNamespace("scTenifoldKnk", quietly = TRUE)) {
  message("scTenifoldKnk ")
  message("devtools::install_github('cailab-tamu/scTenifoldKnk')")
} else {

  library(scTenifoldKnk)
  library(ggrepel)
  set.seed(123)

  knockout_params <- list(
    knockout_genes = "LRRC8D",
    n_hvg = 2000,
    qc_mtThreshold = 0.1,
    qc_minLSize = 1000,
    nc_nNet = 3,
    nc_nCells = 300,
    pval_threshold = 0.05
  )

  message(sprintf("Knockout target gene: %s", paste(knockout_params$knockout_genes, collapse = ", ")))
  message(sprintf("Number of highly variable genes: %d", knockout_params$n_hvg))

  knockout_output_dir <- workDir
  

  scObject_disease_ko <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "Disease")

  message(sprintf("Disease-group cells: %d", ncol(scObject_disease_ko)))
  message(sprintf("Disease-group cell types: %s", paste(unique(Idents(scObject_disease_ko)), collapse = ", ")))

  available_knockout_genes <- knockout_params$knockout_genes[
    knockout_params$knockout_genes %in% rownames(scObject_disease_ko)
  ]
  missing_knockout_genes <- setdiff(knockout_params$knockout_genes, available_knockout_genes)

  if (length(missing_knockout_genes) > 0) {
    message(sprintf("Warning - genes not found: %s", paste(missing_knockout_genes, collapse = ", ")))
  }

  if (length(available_knockout_genes) > 0) {

    message(sprintf(" %d ", length(available_knockout_genes)))

    countMat <- GetAssayData(scObject_disease_ko, layer = "counts")

    scObject_disease_ko <- FindVariableFeatures(
      object = scObject_disease_ko, selection.method = "vst",
      nfeatures = knockout_params$n_hvg
    )
    hvgs <- VariableFeatures(scObject_disease_ko)

    for (target_gene in available_knockout_genes) {

      message(sprintf("\n########## : %s ##########", target_gene))

      gene_output_dir <- workDir
      
      tryCatch({
        data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])
        message(sprintf("  : %d  x %d ", nrow(data), ncol(data)))

        message("  scTenifoldKnk...")
        result <- scTenifoldKnk(
          countMatrix = data,
          gKO = target_gene,
          qc_mtThreshold = knockout_params$qc_mtThreshold,
          qc_minLSize = knockout_params$qc_minLSize,
          nc_nNet = knockout_params$nc_nNet,
          nc_nCells = knockout_params$nc_nCells
        )
        message("  scTenifoldKnk!")

        df <- result$diffRegulation
        df <- df[df$gene != target_gene, ]

        outTab <- df[df$p.adj < knockout_params$pval_threshold, ]
        write.csv(outTab, file = file.path(gene_output_dir, paste0(target_gene, "_sigDiff.csv")),
                  row.names = FALSE)
        message(sprintf("   %d  (p.adj < 0.05)", nrow(outTab)))


        df$log_p.adj <- -log10(df$p.adj)
        df$significant <- ifelse(df$p.adj < knockout_params$pval_threshold, "Significant", "Not significant")
        label_genes <- subset(df, p.adj < knockout_params$pval_threshold)

        y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)

        scatter_colors <- colorRampPalette(brewer.pal(12, "Set3"))(2)
        p_scatter <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
          geom_point(alpha = 0.7, size = 1.5) +
          scale_color_manual(values = c("Significant" = scatter_colors[1], "Not significant" = "gray70")) +
          geom_hline(yintercept = -log10(knockout_params$pval_threshold), linetype = "dashed", color = scatter_colors[1]) +
          geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50,
                          color = "black", fontface = "italic") +
          labs(title = paste0(target_gene, " Knockout"),
               x = "Z-score",
               y = "-log10(p.adj)") +
          theme_classic(base_size = 14) +
          coord_cartesian(ylim = c(0, y_upper)) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 12),
            legend.position = "none"
          )

        pdf(file.path(gene_output_dir, "Figure9G.pdf"), width = 6, height = 5)
        print(p_scatter)
        dev.off()
        message("  : ", target_gene, "_KO_scatter.pdf")


        top_genes <- head(df[order(-df$FC), ], 20)
        bar_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nrow(top_genes))

        p_bar <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC, fill = reorder(gene, FC))) +
          geom_bar(stat = 'identity', alpha = 0.9) +
          scale_fill_manual(values = bar_colors) +
          coord_flip() +
          labs(title = paste0("Top 20 Differentially Regulated Genes\n(", target_gene, " Knockout)"),
               x = "Gene",
               y = "FC") +
          theme_classic(base_size = 14) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(size = 13, face = "bold"),
            axis.text.y = element_text(size = 11, face = "italic"),
            axis.text.x = element_text(size = 11),
            legend.position = "none"
          )

        pdf(file.path(gene_output_dir, "Figure9H.pdf"), width = 6, height = 5)
        print(p_bar)
        dev.off()
        message("  : ", target_gene, "_KO_barplot.pdf")

        write.csv(df, file = file.path(gene_output_dir, paste0(target_gene, "_KO_allResults.csv")),
                  row.names = FALSE)

        message(sprintf("   %s !", target_gene))

      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        message(sprintf("   %s : %s", target_gene, e$message))
      })
    }

  } else {
  }

  message("\n================================================================================")
  message(sprintf("Output directory: %s", knockout_output_dir))

  message("\n===  ===")
  message("      ├── [Gene]_sigDiff.csv         ()")
}

message("\n================================================================================")
message("         ShowGenes Module Score Violin Plot (Disease Group Only)                ")

violin_output_dir <- workDir

tryCatch({
  Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells$cell_type

  valid_genes <- showGenes[showGenes %in% rownames(Peripheral_Blood_Mononuclear_Cells)]
  if (length(valid_genes) == 0) stop("showGenes")
  message(sprintf("Valid hub genes: %d / %d (%s)", length(valid_genes), length(showGenes),
                  paste(valid_genes, collapse = ", ")))

  Peripheral_Blood_Mononuclear_Cells <- AddModuleScore(
    Peripheral_Blood_Mononuclear_Cells,
    features = list(valid_genes),
    name = "ShowGenes_Score"
  )
  message("Module Score")

  disease_obj <- subset(Peripheral_Blood_Mononuclear_Cells, subset = Type == "Disease")
  Idents(disease_obj) <- disease_obj$cell_type
  n_disease <- ncol(disease_obj)
  message(sprintf("Disease-group cells: %d", n_disease))

  unique_cts <- sort(unique(as.character(Idents(disease_obj))))
  n_cts <- length(unique_cts)

  plot_df <- data.frame(
    CellType = as.character(Idents(disease_obj)),
    ModuleScore = disease_obj$ShowGenes_Score1,
    stringsAsFactors = FALSE
  )
  plot_df$CellType <- factor(plot_df$CellType, levels = unique_cts)

  if (n_cts <= 9) {
    ct_cols <- brewer.pal(max(3, n_cts), "Set1")[1:n_cts]
  } else {
    ct_cols <- colorRampPalette(brewer.pal(9, "Set1"))(n_cts)
  }
  names(ct_cols) <- unique_cts

  p_module_vln <- ggplot(plot_df, aes(x = CellType, y = ModuleScore, fill = CellType)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.85, linewidth = 0.4) +
    geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5,
                 outlier.alpha = 0.4, linewidth = 0.4) +
    scale_fill_manual(values = ct_cols) +
    theme_classic(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
    ) +
    labs(y = "Module Score",
         title = "ShowGenes Module Score (Disease Group)",
         subtitle = paste0("Genes: ", paste(valid_genes, collapse = ", ")))

  vln_width <- max(7, n_cts * 0.9 + 2)

  pdf(file.path(violin_output_dir, "Figure9D.pdf"),
      width = vln_width, height = 6)
  print(p_module_vln)
  dev.off()

  write.csv(plot_df, file.path(violin_output_dir, "ModuleScore_Data_Disease.csv"),
            row.names = FALSE)



  indiv_dir <- workDir
  
  expr_disease <- GetAssayData(disease_obj, layer = "data")

  for (g in valid_genes) {
    tryCatch({
      gene_expr <- as.numeric(expr_disease[g, ])
      g_df <- data.frame(
        CellType = factor(as.character(Idents(disease_obj)), levels = unique_cts),
        Expression = gene_expr,
        stringsAsFactors = FALSE
      )

      p_g <- ggplot(g_df, aes(x = CellType, y = Expression, fill = CellType)) +
        geom_violin(trim = FALSE, scale = "width", alpha = 0.85, linewidth = 0.4) +
        geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.5,
                     outlier.alpha = 0.4, linewidth = 0.4) +
        scale_fill_manual(values = ct_cols) +
        theme_classic(base_size = 13) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold.italic", size = 14)
        ) +
        labs(y = "Expression Level", title = g)

      pdf(file.path(indiv_dir, paste0("Supp_", g, "_VlnPlot_Disease.pdf")),
          width = max(7, n_cts * 0.9 + 2), height = 5.5)
      print(p_g)
      dev.off()

    }, error = function(e) {
      if (dev.cur() > 1) dev.off()
      message(sprintf("   %s : %s", g, e$message))
    })
  }

  message(sprintf("Individual plots written to: %s", indiv_dir))

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message(sprintf("ShowGenes Module Score: %s", e$message))
})

message("\n================================================================================")
message("         Ro/e (Observed/Expected Ratio) Heatmap                                 ")


roe_output_dir <- workDir

tryCatch({
  Idents(Peripheral_Blood_Mononuclear_Cells) <- Peripheral_Blood_Mononuclear_Cells$cell_type

  cell_types <- as.character(Idents(Peripheral_Blood_Mononuclear_Cells))
  groups <- Peripheral_Blood_Mononuclear_Cells$Type

  obs_table <- table(cell_types, groups)
  print(obs_table)

  row_sums <- rowSums(obs_table)
  col_sums <- colSums(obs_table)
  total <- sum(obs_table)

  exp_table <- outer(row_sums, col_sums) / total

  roe_table <- obs_table / exp_table

  message("\nRo/e:")
  print(round(roe_table, 3))

  label_table <- matrix("", nrow = nrow(roe_table), ncol = ncol(roe_table))
  for (i in 1:nrow(roe_table)) {
    for (j in 1:ncol(roe_table)) {
      val <- roe_table[i, j]
      if (val > 1.1) {
        label_table[i, j] <- "+"
      } else if (val < 0.9) {
        label_table[i, j] <- "-"
      } else {
        label_table[i, j] <- "±"
      }
    }
  }

  roe_df <- as.data.frame.matrix(roe_table)
  roe_df$CellType <- rownames(roe_df)
  roe_df <- roe_df[, c("CellType", colnames(obs_table))]
  write.csv(roe_df, file.path(roe_output_dir, "Roe_Values.csv"), row.names = FALSE)

  roe_mat <- as.matrix(roe_table)

  roe_long <- data.frame(
    CellType = rep(rownames(roe_mat), ncol(roe_mat)),
    Group = rep(colnames(roe_mat), each = nrow(roe_mat)),
    Roe = as.vector(roe_mat),
    stringsAsFactors = FALSE
  )

  roe_long$Label <- ifelse(roe_long$Roe > 1.1, "+",
                    ifelse(roe_long$Roe < 0.9, "-", "\u00b1"))
  roe_long$DisplayText <- paste0(round(roe_long$Roe, 2), "\n", roe_long$Label)

  roe_long$CellType <- factor(roe_long$CellType, levels = rev(sort(unique(roe_long$CellType))))
  roe_long$Group <- factor(roe_long$Group, levels = c("Control", "Disease"))

  n_ct <- length(unique(roe_long$CellType))
  ht_height <- max(5, n_ct * 0.6 + 2)

  p_roe <- ggplot(roe_long, aes(x = Group, y = CellType, fill = Roe)) +
    geom_tile(color = "grey80", linewidth = 0.5) +
    geom_text(aes(label = DisplayText,
                  color = ifelse(Roe > 2, "white", "black")),
              size = 3.5, lineheight = 0.85) +
    scale_color_identity() +
    scale_fill_gradient2(low = "#F7FBFF", mid = "#9ECAE1", high = "#08519C",
                         midpoint = 1.5, limits = c(0, 3),
                         oob = scales::squish, name = "Ro/e") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
    ) +
    labs(title = "Ro/e: Observed/Expected Cell Ratio (Control vs Disease)")

  pdf(file.path(roe_output_dir, "Figure9F.pdf"),
      width = 6, height = ht_height)
  print(p_roe)
  dev.off()

  message("Ro/e heatmap saved!")

}, error = function(e) {
  if (dev.cur() > 1) dev.off()
  message(sprintf("Ro/e: %s", e$message))
})

sig_file <- file.path(gene_output_dir, paste0(target_gene, "_sigDiff.csv"))
if (file.exists(sig_file)) {
  sig_df <- read.csv(sig_file, header = TRUE, stringsAsFactors = FALSE)
  sig_genes <- unique(sig_df$gene)
  sig_entrez <- tryCatch(
    bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
    error = function(e) character(0))
  if (length(sig_entrez) > 0) {
    kegg_ko <- tryCatch(enrichKEGG(gene = sig_entrez, organism = "hsa",
                                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                        error = function(e) NULL)
    if (!is.null(kegg_ko) && nrow(kegg_ko@result) > 0) {
      ko_res <- kegg_ko@result
      n_show <- min(20, nrow(ko_res))
      p_kegg_dot <- ggplot(ko_res[1:n_show, ],
                           aes(x = Count, y = reorder(Description, Count),
                               color = p.adjust, size = Count)) +
        geom_point() +
        scale_color_gradient(low = "red", high = "blue") +
        labs(title = "KEGG Pathway", x = "Gene Count", y = "") +
        theme_bw(base_size = 12) +
        theme(axis.text.y = element_text(size = 10),
              plot.title = element_text(hjust = 0.5, face = "bold"))
      pdf(file.path(gene_output_dir, "Figure9I.pdf"), width = 9, height = 6)
      print(p_kegg_dot)
      dev.off()
      write.csv(ko_res, file.path(gene_output_dir, paste0(target_gene, "_KO_KEGG.csv")), row.names = FALSE)
    }
  }
}