from scGeneFit.functions import *
import numpy as np
import scanpy as sc
import os
import pandas as pd
np.random.seed(0) 

#Choose an evaluation method (e.g. classification accuracy)
from sklearn.neighbors import NearestCentroid
clf=NearestCentroid()

def performance(X_train, y_train, X_test, y_test, clf):
    clf.fit(X_train, y_train)
    return clf.score(X_test, y_test)

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
parser.add_argument("num_markers", type=int, help="num_markers")
parser.add_argument("eps", type=float, help="eps")
args = parser.parse_args()
    
system = 'spleen'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'

adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system +'.h5ad')
adata = adata[adata.obs["celltype"] != "Unknown", :]

# select HVGs (to minimize complexity)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=10, min_disp=0.2)
adata = adata[:, adata.var.highly_variable]

# read data, celltype labels (as ints) and genes
data = adata.layers["logcounts"].toarray()
genes = adata.var_names
celltypes = adata.obs['celltype'].astype('category').cat.codes.to_numpy()
samples = adata.obs['sample'].astype('category').cat.codes.to_numpy()


# centers
method = 'centers'
current_markers= get_markers(data, celltypes, args.num_markers, method=method, redundancy=0.25, epsilon=args.eps)
current_genes = genes[current_markers]
save_dir = root_dir + 'calibration/spleen_init_settings_scGeneFit/' + method + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
save_file = save_dir + 'nGenes_' + str(len(current_genes)) + '_eps_' + str(current_eps) +'.txt'
with open(save_file, 'w') as f:
    for item in current_genes:
        f.write("%s\n" % item)

# pairwise
method = 'pairwise'
current_markers= get_markers(data, celltypes, args.num_markers, method=method, sampling_rate=0.1, n_neighbors=2, max_constraints=1000, epsilon=args.eps)
current_genes = genes[current_markers]
save_dir = root_dir + 'calibration/spleen_init_settings_scGeneFit/' + method + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
save_file = save_dir + 'nGenes_' + str(len(current_genes)) + '_eps_' + str(current_eps) +'.txt'
with open(save_file, 'w') as f:
    for item in current_genes:
        f.write("%s\n" % item)

# pairwise_centers
method = 'pairwise_centers'
current_markers= get_markers(data, celltypes, args.num_markers, method=method, sampling_rate=0.1, n_neighbors=2, max_constraints=1000, epsilon=args.eps)
current_genes = genes[current_markers]
save_dir = root_dir + 'calibration/spleen_init_settings_scGeneFit/' + method + '/'
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
save_file = save_dir + 'nGenes_' + str(len(current_genes)) + '_eps_' + str(current_eps) +'.txt'
with open(save_file, 'w') as f:
    for item in current_genes:
        f.write("%s\n" % item)  
