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
parser.add_argument("n_samples", type=int, help="n_samples")
parser.add_argument("min_disp", type=float, help="min_disp")
parser.add_argument("n_threads", type=int, help="n_threads")
parser.add_argument("lambda_reg", type=float, help="lambda")
args = parser.parse_args()
    
system = 'spleen'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'
#root_dir = '/Users/alsu/Develop/geneBasis/data/'

adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system + '.h5ad')
samples = adata.obs["sample"].unique()

#n_samples_to_process = list(range(1,(len(samples)+1)))
#min_disp = [0.5, 0.25, 0]
#n_threads = [1,3,6,9]
#lambdas = [2e-4, 4e-4, 6e-4, 8e-4, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3 , 3.5e-3]

save_dir = root_dir + 'processed_time/scmer/' + system + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)

rows = []
t0 = time.time()
current_adata = adata[adata.obs["sample"].isin(samples[0:args.n_samples]),:]
sc.pp.filter_cells(current_adata, min_genes=200)
sc.pp.filter_genes(current_adata, min_cells=3)
sc.pp.normalize_total(current_adata, target_sum=1e4)
sc.pp.log1p(current_adata)
sc.pp.highly_variable_genes(current_adata, min_mean=0.01, max_mean=10, min_disp=args.min_disp)
current_adata = current_adata[:, current_adata.var.highly_variable]
sc.pp.scale(current_adata, max_value=10)
model = scmer.UmapL1(lasso=args.lambda_reg, ridge=0., n_pcs=50, perplexity=100., use_beta_in_Q=True, n_threads=args.n_threads, pca_seed=2020)
model.fit(current_adata.X, batches=current_adata.obs['sample'].values)
current_genes = current_adata.var_names[model.get_mask()].tolist()
t1 = time.time()
t = t1-t0
rows.append([args.n_samples , args.min_disp , args.n_threads , args.lambda_reg , len(current_adata.obs) , len(current_adata.var_names) , len(current_genes), t])
df = pd.DataFrame(rows, columns=["n_processed_samples", "min_disp", "n_threads", "lambda" , "n_cells" , "n_genes_init", "n_genes_selected", "time"])
df.to_csv(save_dir + 'stat_batched_' + str(t1) + '.csv', header=True, index=None, sep=',')
