"""
Converted from multivelo.ipynb to an executable script.
Run with: python multivelo_server.py
"""
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params("scvelo")
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 200)
np.set_printoptions(suppress=True)

# Set all paths for data and output on the cluster
total_counts_path = "/omics/groups/OE0132/tandem/nschmidt/MultiVelo_Data/all_rna_counts_normalized.h5ad"
feature_matrix_data_path = "/omics/groups/OE0132/internal/metzj/projects/Multiome/CellRanger_Output/Exp4/ctr/outs/filtered_feature_bc_matrix/"
peak_annotation = "/omics/groups/OE0132/internal/metzj/projects/Multiome/CellRanger_Output/Exp4/ctr/outs/atac_peak_annotation.tsv"
feature_linkage = "/omics/groups/OE0132/internal/metzj/projects/Multiome/CellRanger_Output/Exp4/ctr/outs/analysis/feature_linkage/feature_linkage.bedpe"

nn_idx_path = "/omics/groups/OE0132/tandem/nschmidt/MultiVelo_Data/nn_idx.txt"
nn_dist_path = "/omics/groups/OE0132/tandem/nschmidt/MultiVelo_Data/nn_dist.txt"
nn_cells_path = "/omics/groups/OE0132/tandem/nschmidt/MultiVelo_Data/nn_cells.txt"

output_path = "/omics/groups/OE0132/tandem/nschmidt/MultiVelo_Data/"

# Begin processing data
adata_rna = sc.read_h5ad(str(total_counts_path))

scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, n_top_genes=1000)
adata_rna.obs.rename(columns={"cell_type": "celltype"}, inplace=True)

adata_atac = sc.read_10x_mtx(feature_matrix_data_path, var_names='gene_symbols', gex_only=False)
adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]


atac_cells_renamed = [cell.replace("-1", "") for cell in adata_atac.obs_names]
adata_atac.obs_names = atac_cells_renamed

rna_cells = adata_rna.obs_names
atac_cells = adata_atac.obs_names

intersection_cells = rna_cells.intersection(atac_cells)

adata_rna_intersected = adata_rna[intersection_cells, :]
adata_atac_intersected = adata_atac[intersection_cells, :]

adata_atac_agg = mv.aggregate_peaks_10x(adata_atac_intersected, str(peak_annotation), str(feature_linkage))

shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac_agg.obs_names))
shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac_agg.var_names))

adata_rna = sc.read_h5ad(str(total_counts_path))
scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, n_top_genes=1000)
adata_rna.obs.rename(columns={"cell_type": "celltype"}, inplace=True)

adata_atac = adata_atac_agg
adata_rna = adata_rna[shared_cells, shared_genes]
adata_atac = adata_atac[shared_cells, shared_genes]

scv.pp.normalize_per_cell(adata_rna)
scv.pp.log1p(adata_rna)
scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)

scv.tl.umap(adata_rna)

nn_idx = np.loadtxt(str(nn_idx_path), delimiter=',')
nn_dist = np.loadtxt(str(nn_dist_path), delimiter=',')
nn_cells = pd.Index(pd.read_csv(str(nn_cells_path), header=None)[0])
nn_cells = nn_cells.str.replace('-1', '', regex=False)

mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)

adata_result = mv.recover_dynamics_chrom(adata_rna,
                                            adata_atac,
                                            max_iter=5,
                                            init_mode="invert",
                                            parallel=True,
                                            save_plot=False,
                                            rna_only=False,
                                            fit=True,
                                            n_anchors=500,
                                            extra_color_key='celltype')

# Save results and compute velocity outputs
adata_result.write(output_path / "multivelo_results.h5ad")
mv.velocity_graph(adata_result)
mv.latent_time(adata_result)
mv.velocity_embedding_stream(adata_result, basis='umap', color='celltype')
plt.savefig(output_path / 'multivelo_velocity_stream_umap.png', dpi=150)
