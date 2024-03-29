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
system = 'atlas_E8.5'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'

# read data
adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system + '.h5ad')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)

# create lambda grid
# lambdas = [1e-4 = 395, 2e-4 = 295, 4e-4 = 180, 6e-4 = 133, 8e-4 = 131, 1e-3 = 116, 1.5e-3 = 84, 2e-3 =45, 2.5e-3=42, 3e-3 , 3.5e-3, 4e-3]
lambdas = [3e-4, 5-4, 1.25e-3, 1.75e-3]
for current_lambda in lambdas:
    model = scmer.UmapL1(lasso=current_lambda, ridge=0., n_pcs=50, perplexity=100., use_beta_in_Q=True, n_threads=9, pca_seed=2020)
    model.fit(adata.X)
    genes = adata.var_names[model.get_mask()].tolist()
    save_dir = root_dir + 'scmer_res/' + system + '/batched/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_file = save_dir + 'nGenes_' + str(len(genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in genes:
            f.write("%s\n" % item)
