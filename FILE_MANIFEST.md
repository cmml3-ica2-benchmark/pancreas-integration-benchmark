# File Manifest

- `ICA2_pancreas_batch_correction_RECOVERED.ipynb`: main documented workflow for loading the processed pancreas dataset, filtering cells, running integrations, benchmarking, and exporting report inputs
- `run_seurat_rpca_preprocessed.R`: Seurat RPCA integration script that reads notebook exports from `seurat_input/` and writes the RPCA embeddings back to CSV
- `make_polished_benchmark_figures_v6.py`: plot-only script for regenerating the final polished figures and supplementary tables from saved `.csv` and `.h5ad` inputs
- `README.md`: practical project overview, software requirements, folder structure, exact run order, and output list
- `requirements.txt`: minimal Python package list for rerunning the notebook and final plotting workflow
- `.gitignore`: excludes large local data files, caches, and other non-essential generated files from a standard GitHub upload
- `results/benchmark_metrics.csv`: source-of-truth core benchmark metrics for the final report
- `results/benchmark_metrics_with_rank.csv`: descriptive ranking table derived from the core benchmark metrics
- `results/Figure2_benchmark_main_wide.csv`: wide-format Figure 2 metrics table derived from `results/benchmark_metrics.csv`
- `results/Figure2_benchmark_main_long.csv`: long-format Figure 2 metrics table derived from `results/benchmark_metrics.csv`
- `results/polished_figure_input.h5ad`: compact AnnData object containing the final UMAP coordinates used by the polished plotting script
- `results/adata_bench_for_scib.h5ad`: benchmark AnnData object containing stored high-dimensional embeddings used for metric calculation and sensitivity analysis
- `results/highdim_resampling/resampling_metrics_summary_highdim_flat.csv`: stored sensitivity summary used to redraw the resampling figure
- `figures_polished_v6/`: final polished figure directory containing the report-ready figures and supplementary outputs
- `seurat_input/`: exported matrix, gene list, barcode list, and metadata used as inputs to the Seurat RPCA script
- `run_seurat_rpca.R`: older RPCA script kept for reference only
- `run_seurat_cca.R`: exploratory Seurat CCA script kept for reference only and not used in the final report

## Article Supplement Mapping

- `figures_polished_v6/Supplementary_Figure_S2_CONCORD_v6.(pdf|png|svg)`: used as Supplementary Figure S1 in the actual article
- `figures_polished_v6/Supplementary_Figure_S3_scib_results_table_scaled_v6.(pdf|png|csv)`: used as Supplementary Figure S2 in the actual article
- `results/highdim_resampling/Supplementary_Figure_S4_resampling_sensitivity_highdim.(pdf|png)`: used as Supplementary Figure S3 in the actual article
- `figures_polished_v6/Supplementary_Table_S1_core_metrics_v6.(csv|pdf|png)`: used as Supplementary Table S1 in the actual article
- `figures_polished_v6/Supplementary_Table_S2A_batch_composition_v6.csv`: part of Supplementary Table S2 in the actual article
- `figures_polished_v6/Supplementary_Table_S2B_celltype_composition_v6.csv`: part of Supplementary Table S2 in the actual article
- `figures_polished_v6/Supplementary_Table_S2C_batch_by_celltype_counts_v6.csv`: part of Supplementary Table S2 in the actual article
