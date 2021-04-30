import scanpy as sc
import os
import pandas as pd
import numpy as np
import pickle as pkl
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats
sc.settings.verbosity = 3
import sys
sys.path.insert(0,'..')
import scmer
import time

system = 'atlas_E8.5'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'
#root_dir = '/Users/alsu/Develop/geneBasis/data/'

adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system + '.h5ad')
samples = adata.obs["sample"].unique()

n_samples_to_process = list(range(1,(len(samples)+1)))
min_disp = [0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5]
n_threads = [1,3,6,9]
lambdas = [1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4, 1e-3, 1.25e-3, 1.5e-3, 1.75e-3, 2e-3, 2.25e-3, 2.5e-3, 2.75e-3, 3e-3 , 3.5e-3, 4e-3]

rows = []
for current_n_samples_to_process in n_samples_to_process:
    for current_min_disp in min_disp:
        for current_n_threads in n_threads:
            for current_lambda in lambdas:
                t0 = time.time()
                current_adata = adata[adata.obs["sample"].isin(samples[0:current_n_samples_to_process]),:]
                sc.pp.filter_cells(current_adata, min_genes=200)
                sc.pp.filter_genes(current_adata, min_cells=3)
                sc.pp.normalize_total(current_adata, target_sum=1e4)
                sc.pp.log1p(current_adata)
                sc.pp.highly_variable_genes(current_adata, min_mean=0.0125, max_mean=3, min_disp=current_min_disp)
                current_adata = current_adata[:, current_adata.var.highly_variable]
                sc.pp.scale(current_adata, max_value=10)
                model = scmer.UmapL1(lasso=current_lambda, ridge=0., n_pcs=50, perplexity=100., use_beta_in_Q=True, n_threads=current_n_threads, pca_seed=2020)
                model.fit(current_adata.X)
                current_genes = current_adata.var_names[model.get_mask()].tolist()
                t1 = time.time()
                t = t1-t0
                rows.append([current_n_samples_to_process , current_min_disp , current_n_threads , current_lambda , len(current_adata.obs) , len(current_adata.var_names) , len(current_genes), t])     
df = pd.DataFrame(rows, columns=["n_processed_samples", "min_disp", "n_threads", "lambda" , "n_cells" , "n_genes_init", "n_genes_selected", "time"])
save_dir = root_dir + 'processed_time/scmer/' + system + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
df.to_csv(save_dir + 'stat_batched.csv', header=True, index=None, sep=',')   
