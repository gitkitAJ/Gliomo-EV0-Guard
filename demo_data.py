"""Generate a small synthetic dataset so the analytics pipeline can be tested
without downloading human genomic data."""

from pathlib import Path
import numpy as np
import pandas as pd

def make_demo(out_dir, seed=42):
    rng = np.random.default_rng(seed)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    patients = [f"P{i:02d}" for i in range(1, 11)]
    rows = []

    genes = ["SOX10", "OLIG2", "NES", "MKI67", "TOP2A", "PDGFRA", "DCX",
             "GFAP", "SLC1A3", "ASCL1", "PROM1", "CDKN1A"]

    for p in patients:
        baseline = rng.normal(0, 0.35)
        for timepoint in ["earlier", "recurrent"]:
            progression = 0 if timepoint == "earlier" else 1
            methylation = 0.75 - 0.18 * progression + baseline + rng.normal(0, .08)
            stem = 0.20 + 0.22 * progression + rng.normal(0, .05)
            differentiated = 0.55 - 0.20 * progression + rng.normal(0, .05)
            cycling = 0.12 + 0.15 * progression + rng.normal(0, .04)

            for cell in range(120):
                rows.append({
                    "patient_id": p,
                    "timepoint": timepoint,
                    "cell_id": f"{p}_{timepoint}_{cell}",
                    "methylation_mean": np.clip(methylation + rng.normal(0,.04), 0, 1),
                    "stem_like_score": stem + rng.normal(0,.07),
                    "differentiated_score": differentiated + rng.normal(0,.07),
                    "cell_cycle_score": cycling + rng.normal(0,.05),
                    "progressed": progression
                })

    cells = pd.DataFrame(rows)

    # Gene-expression-like values for demonstration only.
    for g in genes:
        if g in ["SOX10", "NES", "ASCL1", "PROM1", "MKI67", "TOP2A"]:
            base = cells["progressed"].to_numpy() * 1.5
        elif g in ["GFAP", "SLC1A3", "CDKN1A"]:
            base = (1 - cells["progressed"].to_numpy()) * 1.2
        else:
            base = np.ones(len(cells))
        cells[g] = np.maximum(0, base + rng.normal(0, .5, len(cells)))

    cells.to_csv(out_dir / "demo_cells.csv", index=False)

if __name__ == "__main__":
    make_demo("data/raw")
