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
system = 'melanoma'
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
# lambdas = [1e-4 = 237, 2e-4 = 78, 4e-4 = 15, 6e-4 = 12, 8e-4 = 6, 1e-3 = 0]
lambdas = [1e-5, 5e-5, 1.2e-4, 1.4e-4, 1.6e-4, 1.8e-4, 2.5e-4, 3e-4, 3.5e-4]
for current_lambda in lambdas:
    model = scmer.UmapL1(lasso=current_lambda, ridge=0., n_pcs=50, perplexity=100., use_beta_in_Q=True, n_threads=9, pca_seed=2020)
    model.fit(adata.X, batches=adata.obs['sample'].values)
    genes = adata.var_names[model.get_mask()].tolist()
    save_dir = root_dir + 'scmer_res/' + system + '/unbatched/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_file = save_dir + 'nGenes_' + str(len(genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in genes:
            f.write("%s\n" % item)
