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
import argparse

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("n_samples_discard_blood", type=int, help="n_samples_discard_blood")
parser.add_argument("lambda_reg", type=float, help="lambda")
args = parser.parse_args()
    
system = 'atlas_E8.5'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'
#root_dir = '/Users/alsu/Develop/geneBasis/data/'

adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system + '.h5ad')
samples = adata.obs["sample"].unique()

celltypes_2retract = ["Erythroid1" , "Erythroid2" , "Erythroid3" , "Blood progenitors 1" , "Blood progenitors 2"]

save_dir = root_dir + 'gene_selection_batch_effect/' + system + '/scmer/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

def update_sce(sample, i):
    current_adata = adata[adata.obs["sample"] == sample]
    position_sample = samples_list.index(sample)
    if position_sample < i:
        current_adata = current_adata[~current_adata.obs["celltype"].isin(celltypes_2retract), :]
    return current_adata

# join together

adata_17 = update_sce(17, args.n_samples_discard_blood)
adata_29 = update_sce(29, args.n_samples_discard_blood)
adata_36 = update_sce(36, args.n_samples_discard_blood)
adata_37 = update_sce(37, args.n_samples_discard_blood)

current_adata = adata_17.concatenate(adata_29, adata_36, adata_37)

sc.pp.filter_cells(current_adata, min_genes=200)
sc.pp.filter_genes(current_adata, min_cells=3)
sc.pp.normalize_total(current_adata, target_sum=1e4)
sc.pp.log1p(current_adata)
sc.pp.highly_variable_genes(current_adata, min_mean=0.01, max_mean=10, min_disp=0.2)
current_adata = current_adata[:, current_adata.var.highly_variable]
sc.pp.scale(current_adata, max_value=10)
model = scmer.UmapL1(lasso=args.lambda_reg, ridge=0., n_pcs=50, perplexity=100., use_beta_in_Q=True, n_threads=3, pca_seed=2020)
model.fit(current_adata.X, batches=current_adata.obs['sample'].values)
current_genes = current_adata.var_names[model.get_mask()].tolist()
save_file = save_dir + 'nGenes_' + str(len(current_genes)) + '_n_', args.n_samples_discard_blood , '.txt'
with open(save_file, 'w') as f:
    for item in current_genes:
        f.write("%s\n" % item)
