"""
Download GEO metadata / processed files.

This script intentionally does not download restricted raw sequencing data.
The Nature Genetics paper states that processed data are public and raw data
are access-restricted.

Usage:
    python src/download_geo.py
"""

from pathlib import Path
import GEOparse
import json

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

ACCESSIONS = ["GSE292025", "GSE291885", "GSE292130"]

def main():
    manifest = {}
    for accession in ACCESSIONS:
        print(f"\nDownloading GEO metadata: {accession}")
        gse = GEOparse.get_GEO(geo=accession, destdir=str(RAW), silent=False)

        out = RAW / f"{accession}_metadata.tsv"
        gse.phenotype_data.to_csv(out, sep="\t")

        manifest[accession] = {
            "metadata_file": str(out),
            "n_samples": int(len(gse.phenotype_data)),
            "title": str(gse.metadata.get("title", [""])[0])
        }

    (RAW / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print("\nDone. Inspect data/raw before running the full pipeline.")

if __name__ == "__main__":
    main()
