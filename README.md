# COPD Ion Channel Biomarkers - Analysis Code

R and Python code for identifying ion channel-related biomarkers in chronic
obstructive pulmonary disease (COPD) using integrated bulk and single-cell
transcriptomics, machine learning, and in silico drug screening.

**No data files are distributed with this repository.** All input files must be
downloaded from public databases or prepared by the user as described below.

---

## 1. Repository contents

```
.
├── README.md
├── LICENSE
├── .gitignore
├── Fig2.R  ...  Fig9.R        analysis scripts (run in this order)
├── Fig10.py, Fig10.R          virtual screening and visualization
└── data/
    ├── README.md              input file specifications
    ├── hub_genes_9.csv    nine hub HICGs (provided, see note below)
    ├── methods.csv            107 algorithm combinations (provided)
    ├── ML.R                   machine-learning function library (provided)
    ├── gene set/              HICG gene set (provided)
    └── (all other input files must be added by the user)
```

---

## 2. Data preparation

All input files go into `data/`. They fall into three groups.

### 2.1 Download from public databases

| File | Source | Notes |
|------|--------|-------|
| `data/raw/GSE76925/`<br>`data/raw/GSE20257/`<br>`data/raw/GSE106986/`<br>`data/raw/GSE38974/`<br>`data/raw/GSE8545/` | GEO | Each folder needs the series matrix file (`GSE*_series_matrix.txt`) **and** the platform file (`GPL*.txt`). GSE76925 and GSE20257 form the discovery cohort; GSE106986, GSE38974 and GSE8545 are the three independent validation cohorts. |
| `data/raw_10x/<sample>/` | GEO **GSE173896** | One folder per sample containing the standard 10X files `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`. Folder names must start with `Control` or `Disease` (5 COPD patients, 3 healthy controls). |
| `data/c2.cp.kegg.symbols.gmt` | MSigDB | KEGG gene sets in symbol format (C2:CP). Registration is required; the file is not redistributed here. |
| `data/refer.txt` | CIBERSORT | LM22 signature matrix (tab-delimited, genes in rows, 22 cell types in columns). |

### 2.2 Prepare from the downloaded GEO data

| File | How to build it |
|------|-----------------|
| `data/sample_info.xlsx` | **Created manually by the user** - one file per bulk dataset, placed in the corresponding `data/raw/<GSE>/` folder. It contains two columns with a header row: column 1 = sample accession (e.g. `GSM265019`), column 2 = group label (`con` for control, `tre` for disease). Example:<br>`Accession` \| `Source name`<br>`GSM265019` \| `con`<br>`GSM265020` \| `tre`<br>The group labels are normalised automatically, so `Control`/`Disease` are also accepted. |
| `data/Training_expr.txt`<br>`data/Training_class.txt`<br>`data/Testing_expr.txt`<br>`data/Testing_class.txt` | Machine-learning inputs. `Training_expr.txt` and `Testing_expr.txt`: tab-delimited expression matrices whose first, unnamed column holds the 13 intersection gene symbols and whose header row holds the sample IDs. `Training_class.txt`: two columns, `sample` and `outcome` (0 = control, 1 = COPD), covering 286 samples (the discovery cohort produced by `Fig2.R`). `Testing_class.txt`: three columns, `Sample`, `Cohort` and `outcome`, covering 87 samples from the three validation cohorts (GSE106986 = 19, GSE38974 = 32, GSE8545 = 36). Sample IDs must match the columns of the corresponding expression matrix. |
| `data/test_GSE8545.csv`<br>`data/test_GSE106986.csv`<br>`data/test_GSE38974.csv` | One expression matrix per validation cohort (samples in rows, hub genes in columns, comma-separated; first column `Sample`, IDs ending with `_con` or `_tre`). **Generated automatically into `data/` by `Fig2.R`** - no manual preparation needed. |
| `merged_after_batch_removal.csv` | Batch-corrected expression matrix of the discovery cohort (first column `geneSymbol`, then one column per sample). **Generated automatically by `Fig2.R` in the repository root** (not in `data/`) and reused by `Fig3.R`-`Fig8.R` and `Fig10.py`. Sample names end with `_con` (control) or `_tre` (disease). |
| `data/genes.txt`, `data/probesets.txt` | MCP-counter annotation files. These can be exported from the `MCPcounter` R package, or downloaded from the MCP-counter resource. |
| `data/Inflammatory Cytokines gene.txt` | Plain text, one cytokine gene symbol per line (24 genes). |
| `data/clinical.xlsx` | Clinical variables (first column = sample ID, remaining columns = Gender, Age, Pack-years, etc.). Use de-identified data only. |

### 2.3 Provided in this repository

| File | Description |
|------|-------------|
| `data/hub_genes_9.csv` | The nine hub HICGs, one gene symbol per row with the header `Gene`. |
| `data/methods.csv` | The 107 machine-learning algorithm combinations used in Figure 3, one model per row with the header `Model`. |
| `data/ML.R` | Machine-learning function library (`RunML`, `ExtractVar`, `CalPredictScore`, `PredictClass`, `RunEval`, `scaleData`, `SimpleHeatmap`). |
| `data/gene set/*.txt` | HICG gene set, one column named `Genes`. |

---

## 3. How to run

1. Clone this repository.
2. Copy all required data into `data/` as described above.
3. **Set the R working directory to the repository root** — every script starts with
   `work_dir <- getwd()`, so the scripts must be run from the folder that contains
   `Fig2.R`.
4. Run the scripts in the order below. Output figures and tables are written into
   the repository root.

## 4. Execution order

| Step | Script | Language | Main input | Main output |
|------|--------|----------|-----------|-------------|
| 1 | `Fig2.R` | R | `data/raw/<GSE>/`, `data/gene set/`, `data/sample_info.xlsx` | `merged_after_batch_removal.csv`, `intersection_genes_13.csv`, Figure 2A-H |
| 2 | `Fig3.R` | R | `merged_after_batch_removal.csv`, `data/hub_genes_9.csv`, `data/methods.csv`, `data/ML.R`, `data/Training_*.txt`, `data/Testing_*.txt`, `data/test_GSE*.csv` | Figure 3A-E, `model.rds`, `AUC_mat.txt` |
| 3 | `Fig4.R` | R | `merged_after_batch_removal.csv`, `data/hub_genes_9.csv` | Figure 4A-F, `shap_top_genes.csv` |
| 4 | `Fig5.R` | R | `merged_after_batch_removal.csv`, `data/hub_genes_9.csv`, `data/c2.cp.kegg.symbols.gmt` | Figure 5_<gene>_*.pdf |
| 5 | `Fig6.R` | R | `merged_after_batch_removal.csv`, `data/hub_genes_9.csv`, `data/refer.txt`, `data/genes.txt`, `data/probesets.txt`, `data/Inflammatory Cytokines gene.txt` | Figure 6A-C |
| 6 | `Fig7.R` | R | `merged_after_batch_removal.csv`, `data/hub_genes_9.csv`, `data/clinical.xlsx`, `data/Inflammatory Cytokines gene.txt`, `data/genes.txt`, `data/probesets.txt` | Figure 7A-I, `geneCluster.txt` |
| 7 | `Fig8.R` | R | `merged_after_batch_removal.csv`, `geneCluster.txt` | Figure 8A-B, Figure8_C1/C2_GO_*.pdf, Figure8_C1/C2_KEGG.pdf |
| 8 | `Fig9.R` | R | `data/raw_10x/<sample>/`, `data/hub_genes_9.csv` | Figure S1A-C, Figure 9A-I, `scRNA_matrix.rds` |
| 9 | `Fig10.py` | Python | `merged_after_batch_removal.csv`, DrugReflector | `input.h5ad`, `data/vscore.txt`, `results.csv` |
| 10 | `Fig10.R` | R | `results.csv` | Figure 10A-B |

### Intermediate files passed between scripts

| File | Produced by | Consumed by |
|------|------------|-------------|
| `merged_after_batch_removal.csv` | `Fig2.R` | `Fig3.R` - `Fig8.R`, `Fig10.py` |
| `data/test_GSE*.csv` | `Fig2.R` | `Fig3.R` (panels 3B-3D) |
| `geneCluster.txt` | `Fig7.R` | `Fig8.R` |
| `results.csv` | `Fig10.py` | `Fig10.R` |

---

## 5. Third-party tools

| Tool | Version | Used for |
|------|---------|----------|
| R | 4.5.2 | all statistical analyses |
| Python | 3.13.12 | v-score computation, DrugReflector |
| PyMOL | 2.1.0 | protein preparation, docking visualization |
| AutoDock Tools | 1.5.6 | hydrogen and charge assignment |
| PyRx / AutoDock Vina | 2.0 | molecular docking (blind docking, exhaustiveness = 8) |
| Discovery Studio | 2020 Client | docking interaction analysis |
| DrugReflector | - | drug repurposing prediction (third-party, install separately) |

DrugReflector is not redistributed here. Install it from its official repository,
then run:

```bash
pip install torch numpy pandas scipy anndata scanpy
python predict.py input.h5ad --top-n 500 \
  --model1 checkpoints/model_fold_0.pt \
  --model2 checkpoints/model_fold_1.pt \
  --model3 checkpoints/model_fold_2.pt \
  --output results.csv
```

---

## 6. Statistical analysis

Unless otherwise stated, P values were adjusted for multiple testing using the
Benjamini-Hochberg method, applied **separately within each set of comparisons**;
BH-adjusted P < 0.05 was considered statistically significant.

Differential expression analysis (limma) used a fold-change threshold of
|log2FC| > 0 together with BH-adjusted P < 0.05, the lenient fold-change
threshold being chosen because of the limited number of HICGs.

---

## 7. Notes

### 7.1 Origin of the nine hub genes

`data/hub_genes_9.csv` contains the nine hub HICGs used throughout Figures 3-10.
They were obtained by intersecting the COPD-associated DEGs with the HICG set
(13 genes, listed in `intersection_genes_13.csv`, produced by `Fig2.R`) and then
retaining the genes selected by the machine-learning procedure implemented in
`Fig3.R` (107 algorithm combinations; the genes present in the model with the
best cross-cohort AUC). Users who wish to regenerate this list should run the
107-model pipeline on their own training data and export the selected features.

### 7.2 Molecular docking

The molecular docking step (Figure 10C-K) was carried out **manually** in PyMOL
2.1.0, AutoDock Tools 1.5.6 and PyRx/AutoDock Vina 2.0. No script is provided
for this step because it involves interactive graphical operations. The exact
parameters (receptor models, grid-box centre and size, exhaustiveness, and
binding-energy thresholds) are documented in Supplementary Table S2, and the
resulting docking poses were rendered in PyMOL.

### 7.3 Other notes

- All findings are computational predictions and require experimental validation.
- Figure 3H (co-expression heatmap of the nine hub genes) is produced by a
  separate script that is not included in this repository.
- Figure legends and the manuscript are provided separately.