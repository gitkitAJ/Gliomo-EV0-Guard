import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def save_all(sample_df, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if "timepoint" in sample_df:
        plt.figure(figsize=(9, 6))
        sns.boxplot(data=sample_df, x="timepoint", y="evolution_pressure")
        sns.stripplot(data=sample_df, x="timepoint", y="evolution_pressure",
                      color="black", alpha=.45)
        plt.title("Exploratory evolution-pressure score by timepoint")
        plt.tight_layout()
        plt.savefig(out_dir / "evolution_pressure_timepoint.png", dpi=200)
        plt.close()

    plt.figure(figsize=(9, 6))
    sns.scatterplot(
        data=sample_df,
        x="methylation_loss",
        y="stem_like",
        hue="timepoint" if "timepoint" in sample_df else None
    )
    plt.title("Methylation loss vs stem-like state")
    plt.tight_layout()
    plt.savefig(out_dir / "methylation_vs_stem_like.png", dpi=200)
    plt.close()

    if "patient_id" in sample_df and "timepoint" in sample_df:
        for pid, g in sample_df.groupby("patient_id"):
            if set(["earlier", "recurrent"]).issubset(set(g["timepoint"])):
                plt.plot(g["timepoint"], g["evolution_pressure"], marker="o",
                         alpha=.5)
        plt.title("Patient-level longitudinal trajectories")
        plt.xlabel("Timepoint")
        plt.ylabel("Evolution pressure")
        plt.tight_layout()
        plt.savefig(out_dir / "patient_trajectories.png", dpi=200)
        plt.close()
