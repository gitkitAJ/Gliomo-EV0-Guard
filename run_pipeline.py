import argparse
from pathlib import Path
import pandas as pd

from demo_data import make_demo
from analysis import (
    add_evolution_features,
    sample_summary,
    correlation_table,
    longitudinal_patient_analysis,
    mixed_effects_model
)
from visualize import save_all

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--data_dir", default=str(ROOT / "data" / "raw"))
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    processed = ROOT / "data" / "processed"
    reports = ROOT / "reports"
    processed.mkdir(exist_ok=True)
    reports.mkdir(exist_ok=True)

    if args.demo:
        make_demo(data_dir)

    path = data_dir / "demo_cells.csv"
    if not path.exists():
        raise FileNotFoundError(
            "Expected a harmonized cell table. Run --demo first, or create "
            "data/raw/demo_cells.csv with the required columns."
        )

    cells = pd.read_csv(path)
    cells = add_evolution_features(cells)

    cells.to_csv(processed / "cell_level_features.csv", index=False)

    samples = sample_summary(cells)
    samples.to_csv(processed / "sample_level_summary.csv", index=False)

    corr = correlation_table(cells)
    corr.to_csv(processed / "correlations.csv", index=False)

    longitudinal = longitudinal_patient_analysis(samples)
    longitudinal.to_csv(processed / "longitudinal_test.csv", index=False)

    try:
        model = mixed_effects_model(samples)
        if model is not None:
            (reports / "mixed_effects_summary.txt").write_text(model.summary().as_text())
    except Exception as exc:
        (reports / "mixed_effects_summary.txt").write_text(
            f"Model could not be fit: {exc}"
        )

    save_all(samples, reports)

    print("Pipeline complete.")
    print("See data/processed and reports.")

if __name__ == "__main__":
    main()
