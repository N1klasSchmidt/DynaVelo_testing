import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import ArchVelo as av
import threadpoolctl

print("threadpoolctl.threadpool_info():", threadpoolctl.threadpool_info())

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

data_dir_uploaded = "/omics/groups/OE0132/tandem/nschmidt/ArchVelo_Data/"
 
avel = sc.read_h5ad(f"{data_dir_uploaded}archvelo_final_result.h5ad")
avel.obs['lthsc_root'] = (avel.obs['cell_type'] == 'LT-HSC').astype(float)

av.velocity_graph(avel, vkey="velo_s")

av.latent_time(
    avel,
    vkey='velo_s',        
    root_key='lthsc_root',
    end_key=None                
)

avel.write(f"{data_dir_uploaded}archvelo_final_result_rooted.h5ad")