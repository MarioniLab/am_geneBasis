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
# read all dataset
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'
data = pd.read_csv(root_dir + '4scmer/mouse_embryo/mouse_embryo.csv', header=0)
data = data.T
data.columns = data.iloc[0]
data = data[1:]
data.index.name = "Cell"

# select sample 29
data = data[data["sample"] == 29]
obs = data[['sample', 'celltype']]
data.drop(['sample', 'celltype'], axis=1, inplace=True)
adata = sc.AnnData(data)
for i in obs.columns:
    adata.obs[i] = obs[i]
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


# create nGenes grid
nGenes = [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]
for current_nGenes in nGenes:
    model = scmer.UmapL1.tune(target_n_features = current_nGenes, X=adata.X, ridge=0., n_pcs=None, perplexity=100., use_beta_in_Q=True, n_threads=10, pca_seed=2020)
    genes = adata.var_names[model.get_mask()].tolist()
    save_dir = root_dir + 'scmer_res/mouse_embryo/sample_29/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_file = save_dir + 'nGenes_' + str(len(genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in genes:
            f.write("%s\n" % item)
