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
data = pd.read_csv("/nfs/research1/marioni/alsu/geneBasis/data/4scmer/mouse_embryo/mouse_17.csv", header=0)
data = data.T
data.columns = data.iloc[0]
data = data[1:]
data.index.name = "Cell"
obs = data[['celltype']]
data.drop(['celltype'], axis=1, inplace=True)
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
lambdas = [2.5e-3 , 3e-3]
for current_lambda in lambdas:
    model = scmer.UmapL1(lasso=current_lambda, ridge=0., n_pcs=40, perplexity=100., use_beta_in_Q=True, n_threads=2, pca_seed=2020)
    model.fit(adata.X)
    genes = adata.var_names[model.get_mask()].tolist()
    save_file = '/nfs/research1/marioni/alsu/geneBasis/data/scmer_res/mouse_embryo/mouse_17/nGenes_' + str(len(genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in genes:
            f.write("%s\n" % item)



echo "source /nfs/research1/marioni/alsu/miniconda3/bin/activate scmer && python3 /nfs/research1/marioni/alsu/geneBasis/am_geneBasis/scmer/mouse_embryo/scmer_test.py" | bsub -o /nfs/research1/marioni/alsu/test.out -e /nfs/research1/marioni/alsu/test.err           
