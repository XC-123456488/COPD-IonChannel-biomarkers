import os
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

INPUT_FILE = "merged_after_batch_removal.csv"
OUTPUT_H5AD = "input.h5ad"
OUTPUT_VSCORE = os.path.join("data", "vscore.txt")
DRUGREFLECTOR_OUTPUT = "results.csv"

df = pd.read_csv(INPUT_FILE, index_col=0)
control_samples = [c for c in df.columns if c.endswith("_con")]
treat_samples = [c for c in df.columns if c.endswith("_tre")]
if len(control_samples) == 0 or len(treat_samples) == 0:
    sys.exit("Missing Control or Treat samples")

control_data = df[control_samples].values.T
treat_data = df[treat_samples].values.T
gene_names = df.index.tolist()

mean_control = np.mean(control_data, axis=0)
mean_treat = np.mean(treat_data, axis=0)
var_control = np.var(control_data, axis=0, ddof=0)
var_treat = np.var(treat_data, axis=0, ddof=0)
denominator = np.sqrt(var_control + var_treat) + (var_control + var_treat == 0)
vscores = (mean_treat - mean_control) / denominator

adata_vscores = ad.AnnData(
    X=vscores.reshape(1, -1),
    obs=pd.DataFrame(
        index=["Treat_vs_Control"],
        data={
            "comparison": "Treat_vs_Control",
            "reference": "Control",
            "comparison_group": "Treat",
        },
    ),
    var=pd.DataFrame(index=gene_names, data={"gene_symbol": gene_names}),
)
adata_vscores.write(OUTPUT_H5AD)

results_df = pd.DataFrame({
    "gene": gene_names,
    "v_score_Treat_vs_Control": vscores,
    "mean_control": mean_control,
    "mean_treat": mean_treat,
}).sort_values("v_score_Treat_vs_Control", ascending=False)
results_df.to_csv(OUTPUT_VSCORE, sep="\t", index=False, float_format="%.6f")

command = [
    sys.executable, "predict.py", OUTPUT_H5AD,
    "--top-n", "500",
    "--model1", "checkpoints/model_fold_0.pt",
    "--model2", "checkpoints/model_fold_1.pt",
    "--model3", "checkpoints/model_fold_2.pt",
    "--output", DRUGREFLECTOR_OUTPUT,
]
if os.path.exists("predict.py"):
    subprocess.run(command, check=True)