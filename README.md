# IDH-Glioma Evolution & Treatment-Resistance Analytics

## Project objective

Build an advanced, reproducible data-analytics pipeline around the 2026 Nature Genetics study:

**"Longitudinal changes in DNA methylation in IDH-mutant glioma fuel disease progression through altered cell state differentiation."**

The project is designed to analyze publicly available processed data from the study and quantify relationships among:

- DNA methylation / hypomethylation
- G-CIMP state
- malignant stem-like states
- differentiated states
- cell-cycle activity
- longitudinal progression / recurrence
- cell-state plasticity and transition dynamics
- clinical outcome, where metadata permit

**Important:** this is a research/analytics project, not a clinical decision system. It must not be used to choose treatment for a patient or to claim that a treatment will prevent resistance.

## Research questions

1. Does methylation loss increase with tumor progression/recurrence?
2. Is methylation loss associated with a larger stem-like malignant-cell fraction?
3. Does methylation loss associate with lower differentiation?
4. Which genes/pathways are most strongly associated with the stem-like phenotype?
5. Are PRC2-associated targets enriched among differentially methylated regions?
6. Can a multi-omic "evolution pressure score" distinguish less-progressed from more-progressed samples?
7. Does a longitudinal change in this score precede or accompany changes in cellular state?
8. Can we identify candidate molecular patterns for follow-up experimental validation?

## Study data

The source paper reports 36 tumor samples from 19 patients, including 32 matched samples from 15 patients. The study used bulk WES, snRNA-seq and joint single-nucleus methylation/RNA profiling.

Public processed datasets:

- GSE292025 — snXRBS / DNA methylation
- GSE291885 — Smart-seq2 snRNA-seq
- GSE292130 — 10x Chromium snRNA-seq

The paper states that raw sequencing data are access-restricted, while processed study data are publicly available.

## Project architecture

```text
                PUBLIC GEO DATA
                       |
          +------------+-------------+
          |            |             |
       snRNA-seq    snXRBS        metadata
          |            |             |
          +------------+-------------+
                       |
                 QC / harmonization
                       |
              cell-state annotation
                       |
          +------------+-------------+
          |            |             |
      stem-like    differentiated   cycling
          |            |             |
          +------------+-------------+
                       |
              methylation features
                       |
             longitudinal modeling
                       |
          +------------+-------------+
          |            |             |
     progression    plasticity   resistance
        score         score       hypotheses
          |            |             |
          +------------+-------------+
                       |
                 dashboard/report
```

## Recommended stack

Python:
- pandas
- numpy
- scipy
- scikit-learn
- statsmodels
- scanpy
- anndata
- seaborn
- matplotlib
- plotly
- GEOparse
- requests

Optional:
- gseapy
- networkx
- lifelines
- umap-learn

## Quick start

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
```

Then download/organize data:

```bash
python src/download_geo.py
```

Run the demo pipeline first:

```bash
python src/run_pipeline.py --demo
```

Run against local processed data:

```bash
python src/run_pipeline.py --data_dir data/raw
```

Outputs are written to:

```text
data/processed/
reports/
```

## Core analytical design

### 1. Cell-state scoring

The pipeline computes transparent gene-set scores rather than claiming a clinical label.

Example conceptual gene sets:

- stem-like
- differentiated
- proliferative/cell-cycle

Gene sets should be replaced or expanded using the exact gene lists from the paper/supplementary material when reproducing the publication analysis.

### 2. Methylation score

For each cell/sample:

```text
methylation_loss_score =
    standardized(reference_methylation - observed_methylation)
```

Higher values represent greater relative methylation loss.

The exact paper's methylation classification should be used for publication-level replication.

### 3. Evolution pressure score

A research feature for exploratory analytics:

```text
evolution_pressure =
    z(methylation_loss)
  + z(stem_like_fraction)
  - z(differentiated_fraction)
  + z(cell_cycle_fraction)
```

This is NOT a clinical biomarker. It is an exploratory composite for ranking observations and generating hypotheses.

### 4. Longitudinal analysis

For patients with matched samples:

```text
patient
  -> primary/earlier sample
  -> recurrence/progression sample
```

Use paired tests and mixed-effects models rather than treating every cell as an independent patient observation.

### 5. Resistance hypothesis layer

The project should identify associations such as:

```text
methylation loss
       |
       +--> stem-like state enrichment
       |
       +--> reduced differentiation
       |
       +--> altered state-transition behavior
       |
       +--> possible treatment-resistance mechanism
```

The last arrow is a hypothesis-generating relationship, not proof of causation.

## Reproducibility rules

- Keep patient-level splits intact.
- Never randomly split cells from the same patient across train/test.
- Prefer patient-level bootstrap or leave-one-patient-out validation.
- Correct for multiple testing.
- Report effect sizes and confidence intervals.
- Distinguish correlation from causation.
- Avoid clinical treatment recommendations.
