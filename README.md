![Tests](https://github.com/kandabarau/oncoplotpy/actions/workflows/ci.yml/badge.svg)
![Pylint](https://img.shields.io/badge/pylint-9.58-yellowgreen?logo=python&logoColor=white)
![License](https://img.shields.io/github/license/kandabarau/oncoplotpy)
![Version](https://img.shields.io/badge/version-1.0.0-blue)
# OncoPlotPy

Inspired by [maftools::oncoplot](https://www.bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html#72_Oncoplots), this is a collection of tools to control genomic waterfall layouts in [Matplotlib](https://matplotlib.org/). 

it provides granular control over matrix sorting, marginal bar charts, and clinical annotation tracks.

---

## Workflow

### 1. Data Preparation (`OncoMatrix`)
The `OncoMatrix` class handles the heavy lifting of pivoting raw MAF data and managing hierarchical sorting.

* **`prepare()`**: Pivots data into a Gene x Sample matrix.
* **`sort_genes()`**: Ranks genes by mutation frequency.
* **`sort_samples()`**: Executes the waterfall logic. Use the `sortby` parameter to prioritize clinical annotations (e.g., Subtype) before the mutation staircase.



### 2. Plotting Suite
The library splits the oncoplot into modular components that can be used independently or through the `run_oncoplot` wrapper:

* **`plot_waterfall`**: The central mutation heatmap.
* **`plot_tmb`**: Top bar plot for sample mutation burden.
* **`plot_genes`**: Right-side bar plot for gene mutation frequency.
* **`plot_annot`**: Bottom tracks for clinical/categorical data.

---

## Layout Alignment

To maintain perfect alignment across the subplot grid, the internal matrices follow this structure:

| Component | Rows | Columns | Alignment |
| :--- | :--- | :--- | :--- |
| **TMB** | Mutation Types | **Samples** | Heatmap X-axis |
| **Waterfall** | **Genes** | **Samples** | Central Grid |
| **Gene Summary** | **Genes** | Mutation Types | Heatmap Y-axis |



---

## Quick Start

```python
from oncoplotpy import OncoMatrix, run_oncoplot

# Prepare data to follow or adjust default logic
om = OncoMatrix(maf_df, annotation_columns)
om.prepare().sort_genes().sort_samples(sortby=['Subtype'])

# Or simply render in one go
run_oncoplot(om, annot_cols=['Subtype'])
