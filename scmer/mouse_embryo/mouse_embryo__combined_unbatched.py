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

# split by samples
samples = data['sample'].unique()
# create nGenes grid
nGenes = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]
obs = current_data[['sample', 'celltype']]
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
for current_nGenes in nGenes:
    model = scmer.UmapL1.tune(target_n_features=current_nGenes, n_pcs=40, perplexity=100., X=adata.X)
    selected_adata = model.transform(adata)
    genes = selected_adata.var_names.to_list()
    save_dir = root_dir + 'scmer_res/mouse_embryo/combined_unbatched/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_file = save_dir + 'nGenes_' + str(len(genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in genes:
            f.write("%s\n" % item)
