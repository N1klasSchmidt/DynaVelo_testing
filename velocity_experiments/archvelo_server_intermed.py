import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import seaborn as sns
import ArchVelo as av

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

data_dir_uploaded = "/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/"
model_outdir = f"{data_dir_uploaded}archvelo/modeling_results/"

n_neigh = 50
n_pcs = 50
num_comps = 10
n_jobs = 8

# Load processed RNA
adata_rna = sc.read_h5ad("/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/adata_rna_archvelo_proc.h5ad")
# Denoised ATAC (not needed when MultiVelo is already run)
#atac_AA = sc.read_h5ad("/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/archvelo/modeling_results/adata_atac_AA_denoised.h5ad")
# Archetypes
smooth_arch = sc.read_h5ad("/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/archvelo/modeling_results/arches.h5ad")
# Archetypal weights per gene
gene_weights = pd.read_csv("/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/archvelo/modeling_results/gene_weights.csv", index_col = [0])
# MultiVelo output
full_mv_res_denoised = sc.read_h5ad("/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/archvelo/modeling_results/multivelo_result_denoised_chrom.h5ad")

def apply_ArchVelo_intermed(adata_rna, 
                        #atac_AA_denoised,
                        full_res_denoised,
                        smooth_arch, 
                        gene_weights, 
                        model_outdir,
                        gene_list=None, 
                        method='Nelder-Mead', 
                        maxiter1=1500, 
                        max_outer_iter=3, 
                        update_mode='cells', 
                        n_jobs=-1, 
                        n_neighbors=50, 
                        n_pcs=50, 
                        verbose=False):
    
    avel = av.apply_ArchVelo(adata_rna, full_res_denoised, smooth_arch, gene_weights, model_outdir, gene_list=gene_list, method=method, maxiter1=maxiter1, max_outer_iter=max_outer_iter, update_mode=update_mode, n_jobs=n_jobs, verbose=verbose)
    return avel


# Main ArchVelo method, using all the pre-computed results we already have!
avel = apply_ArchVelo_intermed(adata_rna,
                    full_mv_res_denoised,
                    smooth_arch,
                    gene_weights,
                    model_outdir,
                    n_jobs = n_jobs,
                    n_neighbors = n_neigh,
                    n_pcs = n_pcs)

avel.write(f"{data_dir_uploaded}archvelo_intermed_result.h5ad")

del adata_rna
#del atac_AA
del full_mv_res_denoised
del smooth_arch
del gene_weights

av.velocity_graph(avel)
av.latent_time(avel)
# av.velocity_embedding_stream(avel, 
#                              show=False, 
#                              color = 'cell_type_abbr', 
#                              title = 'ArchVelo result')

avel.write(f"{data_dir_uploaded}archvelo_final_result.h5ad")