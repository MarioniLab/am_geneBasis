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

method = 'pairwise_centers'
sampling_rate=0.1
n_neighbors=0
max_constraints=1000

# read dataset
system = 'pancreas'
root_dir = '/nfs/research1/marioni/alsu/geneBasis/data/'
adata = sc.read_h5ad(root_dir + 'scRNA_datasets/' + system + '/sce_' + system +'.h5ad')

# select HVGs (to minimize complexity)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.01, max_mean=3, min_disp=0.01)
adata = adata[:, adata.var.highly_variable]

# read data, celltype labels (as ints) and genes
data = adata.layers["logcounts"].toarray()
genes = adata.var_names
celltypes = adata.obs['celltype'].astype('category').cat.codes.to_numpy()
samples = adata.obs['sample'].astype('category').cat.codes.to_numpy()

# split into train and test to optimize eps
idx_test = np.asarray(np.where(samples == 0))
idx_train = np.asarray(np.where(samples != 0))

data_train = data[idx_train,:]
data_train = data_train[0,:,:]
celltypes_train = celltypes[idx_train].reshape(-1)

data_test = data[idx_test,:]
data_test = data_test[0,:,:]
celltypes_test = celltypes[idx_test].reshape(-1)

# run
num_markers = np.arange(10,260,10)
for current_num_markers in num_markers:
    res=optimize_epsilon(data_train, celltypes_train, data_test, celltypes_test, current_num_markers, method=method, fixed_parameters={'redundancy': 0.1, 'sampling_rate':0.1, 'max_constraints': 200}, 
                        bounds=[(0.1 , 5)], x0=[1], max_fun_evaluations=25, n_experiments=10, 
                        clf=clf, verbose=True)
    current_markers= get_markers(data, celltypes, current_num_markers, method=method, sampling_rate=sampling_rate,
                     n_neighbors=n_neighbors, epsilon=float(res[0]), max_constraints=max_constraints)
    current_genes = genes[current_markers]
    save_dir = root_dir + 'scGeneFit_res/' + system + '/' + method + '/'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    save_file = save_dir + 'nGenes_' + str(len(current_genes)) + '.txt'
    with open(save_file, 'w') as f:
        for item in current_genes:
            f.write("%s\n" % item)
