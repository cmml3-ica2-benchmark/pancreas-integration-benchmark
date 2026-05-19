# Pancreas Batch-Correction Mini-Project

This project benchmarks batch-correction methods on a processed human pancreas scRNA-seq dataset. The final comparison focuses on `Uncorrected`, `Harmony`, `Scanorama`, `Seurat RPCA`, and exploratory `CONCORD`, using batch ASW, batch mixing score, cell-type ASW, Leiden ARI, and Leiden NMI.

The project is organized so the fixed report outputs can be reproduced without changing the scientific analysis logic. `results/benchmark_metrics.csv` is the source of truth for final core metric values when other derived files disagree.

## Requirements
- Python 3.10+ with `scanpy`, `anndata`, `omicverse`, `matplotlib`, `numpy`, `pandas`, `scikit-learn`, and `openpyxl`
- R with `Seurat` and `Matrix`
- Optional Python package for exploratory analysis: `concord`
- A minimal reproducibility file is provided as `requirements.txt`
- Environment snapshots for recovery are stored in `env_recovery/`

## GitHub Repository Contents
- Recommended upload set: `README.md`, `FILE_MANIFEST.md`, `ICA2_pancreas_batch_correction_RECOVERED.ipynb`, `run_seurat_rpca_preprocessed.R`, `make_polished_benchmark_figures_v6.py`, `requirements.txt`, `results/benchmark_metrics.csv`, `results/Figure2_benchmark_main_wide.csv`, `results/Figure2_benchmark_main_long.csv`, and selected final figures if desired
- `results/polished_figure_input.h5ad` is about 270 MB and `results/adata_bench_for_scib.h5ad` is about 281 MB, so they should not be pushed to a standard GitHub repository without Git LFS
- If the `.h5ad` files are kept local, document that the workflow still runs in the order listed below once those files are placed back into `results/`
- The processed input `data/pancreas.h5ad` may also be omitted from GitHub if redistribution is not appropriate; in that case keep the file path and data note in this README
- This repository focuses on code, workflow documentation, and source metric tables; final figures used in the report were curated locally and are not required for the GitHub upload.

## Expected Structure
- `data/pancreas.h5ad`: processed pancreas AnnData used as the notebook input
- `seurat_input/`: matrix and metadata exported from the notebook for Seurat RPCA
- `results/`: source metrics, figure-ready CSVs, compact AnnData objects, and stored sensitivity outputs
- `figures/`: earlier notebook-generated figures
- `figures_polished_v6/`: final polished report figures and supplementary tables
- `results/highdim_resampling/`: stored resampling summaries and sensitivity figure outputs

## Exact Run Order
1. Open `ICA2_pancreas_batch_correction_RECOVERED.ipynb` and run it in order.
2. Stop at the Seurat export step and run `run_seurat_rpca_preprocessed.R`.
3. Return to the notebook and import the Seurat RPCA outputs.
4. Continue the notebook through metric calculation, stored-embedding sensitivity analysis, and final output export.
5. Run `make_polished_benchmark_figures_v6.py` if you need to regenerate the polished figure set from saved inputs only.

## Main Files
- `ICA2_pancreas_batch_correction_RECOVERED.ipynb`: documented end-to-end workflow for loading data, filtering, baseline construction, Harmony, Scanorama, Seurat RPCA import, exploratory CONCORD, benchmarking, sensitivity analysis, and final output export
- `run_seurat_rpca_preprocessed.R`: R script that reads `seurat_input/`, runs Seurat RPCA, and exports `seurat_rpca_pca.csv`, `seurat_rpca_umap.csv`, and `seurat_rpca_metadata.csv`
- `make_polished_benchmark_figures_v6.py`: plot-only script that redraws the final polished figures and tables from saved `.h5ad` and `.csv` files
- `results/benchmark_metrics.csv`: source-of-truth core benchmark metrics used by the report
- `results/benchmark_metrics_with_rank.csv`: descriptive ranking table derived from the same benchmark metrics
- `results/polished_figure_input.h5ad`: compact AnnData object with the final UMAP coordinates used by the polished figure script
- `results/Figure2_benchmark_main_wide.csv`: Figure 2 summary table derived from `results/benchmark_metrics.csv`
- `results/Figure2_benchmark_main_long.csv`: long-format Figure 2 summary table derived from `results/benchmark_metrics.csv`
- `results/highdim_resampling/resampling_metrics_summary_highdim_flat.csv`: stored summary used for the sensitivity figure

## Final Outputs
- `figures_polished_v6/Figure_1_visual_atlas_v6.(png|pdf|svg)`
- `figures_polished_v6/Figure_2_benchmark_synthesis_v6.(png|pdf|svg)`
- `figures_polished_v6/Supplementary_Figure_S1_full_celltype_atlas_v6.(png|pdf|svg)`
- `figures_polished_v6/Supplementary_Figure_S2_CONCORD_v6.(png|pdf|svg)`
- `figures_polished_v6/Supplementary_Figure_S3_scib_results_table_scaled_v6.(png|pdf|csv)`
- `results/highdim_resampling/Supplementary_Figure_S4_resampling_sensitivity_highdim.(png|pdf)`
- `figures_polished_v6/Supplementary_Table_S1_core_metrics_v6.(csv|png|pdf)`
- `figures_polished_v6/Supplementary_Table_S2_dataset_composition_v6.(png|pdf)`
- `results/Supplementary_Table_S1_full_benchmark.(csv|xlsx|png)`

## Article Mapping
- Main Figure 1 in the project files is `figures_polished_v6/Figure_1_visual_atlas_v6.*`
- Main Figure 2 in the project files is `figures_polished_v6/Figure_2_benchmark_synthesis_v6.*`
- Article Supplementary Figure S1 corresponds to `figures_polished_v6/Supplementary_Figure_S2_CONCORD_v6.*`
- Article Supplementary Figure S2 corresponds to `figures_polished_v6/Supplementary_Figure_S3_scib_results_table_scaled_v6.*`
- Article Supplementary Figure S3 corresponds to `results/highdim_resampling/Supplementary_Figure_S4_resampling_sensitivity_highdim.*`
- Article Supplementary Table S1 corresponds to `figures_polished_v6/Supplementary_Table_S1_core_metrics_v6.(csv|png|pdf)`
- Article Supplementary Table S2 is composed of:
- `figures_polished_v6/Supplementary_Table_S2A_batch_composition_v6.csv`
- `figures_polished_v6/Supplementary_Table_S2B_celltype_composition_v6.csv`
- `figures_polished_v6/Supplementary_Table_S2C_batch_by_celltype_counts_v6.csv`

## Notes
- `CONCORD` is exploratory because it uses a 100-dimensional cosine representation rather than the 20-dimensional Euclidean setting used for the main comparison.
- The sensitivity analysis uses stored embeddings and does not represent full reintegration robustness.
- `run_seurat_cca.R` and the Seurat CCA block are kept only for reference and are not part of the final report workflow.
- Older historical figure folders and older plotting scripts may contain pre-v6 intermediate values; use the v6 files and `results/benchmark_metrics.csv` for the final report.
- For GitHub submission, `results/benchmark_metrics.csv` and the derived Figure 2 CSVs are small enough to upload directly and should be preferred over screenshots or copied numbers.
