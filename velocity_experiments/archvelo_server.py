import os
import scipy
import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc
import scvelo as scv
import multivelo as mv

import matplotlib.pyplot as plt
import seaborn as sns

import ArchVelo as av

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

# Load data and set params
data_dir = "/omics/groups/OE0132/internal/metzj/projects/Multiome/CellRanger_Output/Exp4/ctr/outs/"
data_dir_uploaded = "/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/"
adata_rna = sc.read_h5ad(f"{data_dir_uploaded}all_rna_counts.h5ad")
adata_atac_raw = sc.read_h5ad(f"{data_dir_uploaded}all_atac_peaks.h5ad")

peak_annotation_path = rf"{data_dir}atac_peak_annotation.tsv"
peak_annotation = pd.read_csv(peak_annotation_path, sep = '\t')
feature_linkage_path = f"{data_dir}analaysis/feature_linkage/feature_linkage.bedpe"

n_neigh = 50
n_pcs = 50

# for ArchVelo
data_outdir = f"{data_dir_uploaded}archvelo/processed_data/"
model_outdir = f"{data_dir_uploaded}archvelo/modeling_results/"
num_comps = 10
n_jobs = 100

# Prepare ATAC
adata_atac_raw.layers['raw_counts'] = adata_atac_raw.X
sc.pp.filter_cells(adata_atac_raw, min_counts=3000)
sc.pp.filter_genes(adata_atac_raw, min_cells=0.01*adata_atac_raw.shape[0])
raw = adata_atac_raw.layers['raw_counts'].copy()
raw.data = np.ceil(raw.data / 2)
adata_atac_raw.layers['poisson_corrected'] = raw
adata_atac_raw.X = adata_atac_raw.layers["poisson_corrected"].toarray()
sc.experimental.pp.normalize_pearson_residuals(adata_atac_raw, theta=1)
adata_atac_raw.layers["pearson"] = adata_atac_raw.X.copy()
adata_atac_raw_proc = adata_atac_raw.copy()

# Prepare RNA
scv.pp.filter_and_normalize(adata_rna, min_cells_u = 20, min_cells = 20, min_shared_counts = 100)
sc.pp.highly_variable_genes(adata_rna, n_top_genes = 1500, 
                            subset = True, flavor = 'seurat')
sc.pp.scale(adata_rna)

sc.tl.pca(adata_rna, n_pcs)   
sc.pp.neighbors(adata_rna, n_neigh, n_pcs = n_pcs)
sc.tl.umap(adata_rna)

# Subset peaks
chromosome = peak_annotation.loc[:,"chrom"].values
start = peak_annotation.loc[:,"start"].values
end = peak_annotation.loc[:,"end"].values

combined_peak_names = [f"{chromosome[i]}:{start[i]}-{end[i]}" for i in range(len(chromosome))]
peak_annotation["peak_name"] = combined_peak_names
peak_annotation.set_index("peak_name", inplace = True)

peak_annotation = peak_annotation.loc[adata_atac_raw_proc.var_names,:]

# Subset the peaks, for which the corresponding gene is present in the RNA data
mask = peak_annotation["gene"].isin(adata_rna.var_names.to_list())
rel_peaks = peak_annotation[mask].index.values
adata_atac_raw_proc = adata_atac_raw_proc[:, rel_peaks].copy()

# Subset the UNIQUE genes, which are used to map the peaks to the RNA data
rel_genes = np.unique(peak_annotation[peak_annotation['gene'].isin(adata_rna.var_names.to_list()).values]['gene'])

# Multivelo preprocessing
adata_atac_raw_multi = sc.read_h5ad(fr"{data_dir_uploaded}all_atac_peaks.h5ad")

adata_atac_agg_peaks = mv.aggregate_peaks_10x(adata_atac_raw_multi, 
                                    peak_annotation_path, 
                                    feature_linkage_path, 
                                    verbose=True)

mv.tfidf_norm(adata_atac_agg_peaks)

# Compute shared genes
adata_atac_agg_peaks.obs_names = [i[:-2] for i in adata_atac_agg_peaks.obs_names]
adata_atac_raw_proc.obs_names = [i[:-2] for i in adata_atac_raw_proc.obs_names]
shared_cells = pd.Index(np.intersect1d(np.intersect1d(adata_rna.obs_names, adata_atac_agg_peaks.obs_names), adata_atac_raw_proc.obs_names))
shared_genes_atac = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac_agg_peaks.var_names))
shared_genes_atac_raw = pd.Index(np.intersect1d(adata_rna.var_names, rel_genes))

#mapping somewhere
mapped_genes = np.unique(np.unique(peak_annotation['gene'].dropna()))
#present in rna
rna_genes = adata_rna.var_names
#present in adata_atac
aggr_genes = adata_atac_agg_peaks.var_names
#present in adata_atac_raw
atac_raw_genes = np.unique(peak_annotation.loc[adata_atac_raw_proc.var_names]['gene'])

shared_genes_total = mapped_genes
shared_genes_total = np.intersect1d(shared_genes_total, rna_genes)
shared_genes_total = np.intersect1d(shared_genes_total, aggr_genes)
shared_genes_total = np.intersect1d(shared_genes_total, atac_raw_genes)

peaks_and_genes_in_adata = peak_annotation.loc[adata_atac_raw_proc.var_names, "gene"]

peak_mask = peaks_and_genes_in_adata.isin(shared_genes_total).values  # len 32186, of which 24492 True
peaks_in_shared_adata = peaks_and_genes_in_adata.loc[peak_mask].index
rel_peaks_total = adata_atac_raw_proc.var_names.intersection(peaks_in_shared_adata)
adata_atac_raw_proc.var_names_make_unique()
adata_atac_agg_peaks.var_names_make_unique()
adata_rna_intersect = adata_rna[shared_cells, shared_genes_total]
adata_atac_agg_intersect = adata_atac_agg_peaks[shared_cells, shared_genes_total]
adata_atac_raw_proc_intersect = adata_atac_raw_proc[shared_cells, rel_peaks_total]

# Further Process RNA
sc.pp.pca(adata_rna_intersect, n_pcs)   
sc.pp.neighbors(adata_rna_intersect, n_neigh, n_pcs = n_pcs)   
sc.tl.umap(adata_rna_intersect, n_components = 2)
scv.pp.moments(adata_rna_intersect, n_pcs=n_pcs, n_neighbors=n_neigh)

np.random.seed(57)
pal = list(np.array(sns.color_palette('husl', 11))[np.random.choice(11,11, replace = False)])
sc.pl.umap(adata_rna_intersect, color = ['cell_type'], palette = pal)

# Run ArchVelo

# Collect all necessary input data for modeling
adata_rna = adata_rna_intersect.copy()
adata_atac_raw = adata_atac_raw_proc_intersect.copy()

peak_annotation.to_csv(rf"{data_dir_uploaded}atac_peak_annotation_archvelo_proc.tsv", sep = '\t')
peak_annotation_proc = pd.read_csv(rf"{data_dir_uploaded}atac_peak_annotation_archvelo_proc.tsv", sep = '\t')

# Run modeling pipeline
XC_raw, S_raw = av.apply_AA_no_test(adata_atac_raw, k = num_comps,
                  outdir = f"{model_outdir}archetypes/")

_, gene_weights = av.annotate_and_summarize(S_raw, peak_annotation_proc, 
                                            outdir = model_outdir)

gene_weights = gene_weights.loc[:, adata_rna.var_names]

atac_AA = av.create_denoised_atac(adata_rna, gene_weights, 
                                  XC_raw, model_outdir = model_outdir,
                                  n_pcs=n_pcs, n_neighbors=n_neigh)

smooth_arch = sc.read_h5ad(f"{model_outdir}arches.h5ad")

# Main ArchVelo method
avel = av.apply_ArchVelo_full(adata_rna,
                    atac_AA,
                    smooth_arch,
                    gene_weights,
                    model_outdir,
                    n_jobs = n_jobs,
                    n_neighbors = n_neigh,
                    n_pcs = n_pcs)

av.velocity_graph(avel)
av.latent_time(avel)
# av.velocity_embedding_stream(avel, 
#                              show=False, 
#                              color = 'cell_type_abbr', 
#                              title = 'ArchVelo result')

avel.write(f"{data_dir_uploaded}archvelo_result.h5ad")