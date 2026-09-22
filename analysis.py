"""
Core analytics for IDH-glioma evolution.

The functions operate on a cell-level table. For real study data, create a
harmonized table with at least:

patient_id, sample_id, timepoint, cell_id,
methylation_mean,
stem_like_score,
differentiated_score,
cell_cycle_score

Optional:
tumor_type, grade, G_CIMP, treatment, recurrence, outcome
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from sklearn.preprocessing import StandardScaler
import statsmodels.formula.api as smf

def zscore(s):
    s = pd.Series(s)
    return (s - s.mean()) / (s.std(ddof=0) + 1e-9)

def add_evolution_features(df):
    df = df.copy()

    # Relative methylation-loss feature.
    df["methylation_loss_score"] = -zscore(df["methylation_mean"])

    df["evolution_pressure_score"] = (
        zscore(df["methylation_loss_score"])
        + zscore(df["stem_like_score"])
        - zscore(df["differentiated_score"])
        + zscore(df["cell_cycle_score"])
    )

    return df

def sample_summary(df):
    group_cols = [c for c in ["patient_id", "sample_id", "timepoint",
                              "tumor_type", "grade", "G_CIMP", "progressed"]
                  if c in df.columns]

    if not group_cols:
        group_cols = ["patient_id", "timepoint"]

    return (df.groupby(group_cols, dropna=False)
              .agg(
                  n_cells=("cell_id", "nunique"),
                  methylation_mean=("methylation_mean", "mean"),
                  methylation_loss=("methylation_loss_score", "mean"),
                  stem_like=("stem_like_score", "mean"),
                  differentiated=("differentiated_score", "mean"),
                  cell_cycle=("cell_cycle_score", "mean"),
                  evolution_pressure=("evolution_pressure_score", "mean")
              )
              .reset_index())

def correlation_table(df):
    metrics = [
        "methylation_loss_score",
        "stem_like_score",
        "differentiated_score",
        "cell_cycle_score",
        "evolution_pressure_score"
    ]

    out = []
    for i, a in enumerate(metrics):
        for b in metrics[i+1:]:
            x = df[a].dropna()
            y = df.loc[x.index, b].dropna()
            common = x.index.intersection(y.index)
            if len(common) < 10:
                continue
            rho, p = spearmanr(df.loc[common, a], df.loc[common, b])
            out.append({"feature_a": a, "feature_b": b,
                        "spearman_rho": rho, "p_value": p,
                        "n": len(common)})
    return pd.DataFrame(out)

def longitudinal_patient_analysis(sample_df):
    """Paired analysis when each patient has earlier and recurrent samples."""
    required = {"patient_id", "timepoint", "evolution_pressure"}
    if not required.issubset(sample_df.columns):
        return pd.DataFrame()

    pivot = sample_df.pivot_table(
        index="patient_id",
        columns="timepoint",
        values="evolution_pressure",
        aggfunc="mean"
    )

    if not {"earlier", "recurrent"}.issubset(pivot.columns):
        return pd.DataFrame()

    paired = pivot[["earlier", "recurrent"]].dropna()
    if len(paired) < 3:
        return pd.DataFrame()

    stat, p = wilcoxon(paired["earlier"], paired["recurrent"])
    return pd.DataFrame([{
        "n_patients": len(paired),
        "wilcoxon_statistic": stat,
        "p_value": p,
        "median_change": (paired["recurrent"] - paired["earlier"]).median()
    }])

def mixed_effects_model(sample_df):
    """Patient-aware model. Exploratory only."""
    needed = {"patient_id", "evolution_pressure", "progressed"}
    if not needed.issubset(sample_df.columns):
        return None

    model = smf.mixedlm(
        "evolution_pressure ~ progressed",
        sample_df,
        groups=sample_df["patient_id"]
    )
    return model.fit(reml=False)
