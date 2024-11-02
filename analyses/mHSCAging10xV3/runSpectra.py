import h5py
import re
import copy
import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import seaborn as sns
from tqdm.auto import trange, tqdm
import matplotlib.pyplot as plt
import sys
import os
from spectra import spectra_util as spc_tl
from spectra import spectra as spc

def check_gene_set_dictionary(adata, annotations, obs_key='cell_type',global_key='global', return_dict = True):
    '''
    Filters annotations dictionary contains only genes contained in the adata. 
    Checks that annotations dictionary cell type keys and adata cell types are identical.
    Checks that all gene sets in annotations dictionary contain >2 genes after filtering.
    
    adata: AnnData , data to use with Spectra
    annotations: dict , gene set annotations dictionary to use with Spectra
    obs_key: str , column name for cell type annotations in adata.obs
    global_key: str , key for global gene sests in gene set annotation dictionary
    return_dict: bool , return filtered gene set annotation dictionary
    
    returns: dict , filtered gene set annotation dictionary
    
    '''
    #test if keys match
    adata_labels  = list(set(adata.obs[obs_key]))+['global']#cell type labels in adata object
    annotation_labels = list(annotations.keys())
    matching_celltype_labels = list(set(adata_labels).intersection(annotation_labels))
    if set(annotation_labels)==set(adata_labels):
        print('Cell type labels in gene set annotation dictionary and AnnData object are identical')
        dict_keys_OK = True
    if len(annotation_labels)<len(adata_labels):
        print('The following labels are missing in the gene set annotation dictionary:',set(adata_labels)-set(annotation_labels))
        dict_keys_OK = False
    if len(adata_labels)<len(annotation_labels):
        print('The following labels are missing in the AnnData object:',set(annotation_labels)-set(adata_labels))
        dict_keys_OK = False
        
    #check that gene sets in dictionary have len >2
    Counter = 0
    annotations_new = {}
    for k,v in tqdm(annotations.items()):
        annotations_new[k] = {}
        for k2,v2 in v.items():
            annotations_new[k][k2]= [x for x in v2 if x in adata.var_names]
            length = len(v2)
            if length<3:
                print('gene set',k2,'for cell type',k,'is of length',length)
                Counter = Counter+1
            
    if Counter > 0:
        print(Counter,'gene sets are too small. Gene sets must contain at least 3 genes')
    elif Counter == 0 and dict_keys_OK:
        print('Your gene set annotation dictionary is correctly formatted.')
    if return_dict:
        return annotations_new
    
# Load scRNA data
print("Loading scRNA data")
h5_path = "../../data/mHSCAging10xMultiome/scRNA.h5"
hf = h5py.File(h5_path, 'r')
RNA_matrix = hf.get('RNA_matrix')[:, :]
UMAP = hf.get('UMAP')[:, :]
barcodes = hf.get('barcodes')
genes = hf.get('genes')
barcodes = [bc.decode('ascii') for bc in barcodes]
genes = [g.decode('ascii') for g in genes]
hf.close()

# Initialize AnnData object
print("Initializing AnnData object")
adata = ad.AnnData(RNA_matrix)
adata.obs_names = barcodes
adata.var_names = genes
adata.obs["cell_type"] = ["HSC"] * len(barcodes)
adata.obsm["X_umap"] = UMAP.T
adata.obsm["ATAC_umap"] = UMAP.T

# Get age labels
age = [re.split("-", i)[1] for i in list(adata.obs.index)]
adata.obs['age'] = age

# Load annotations
annotations = spc_tl.get_default_dict()
annotations["HSC"] = {}
marker_dir = "../../data/mHSCAging10xMultiome/markers/RNA/"
marker_files = [i for i in os.listdir(marker_dir) if re.search("txt", i)]
marker_names = [i.replace(".txt", "") for i in marker_files]
for marker_ind in range(len(marker_names)):
    annotations["HSC"][marker_names[marker_ind]] = list(pd.read_table(
        "../../data/mHSCAging10xMultiome/markers/RNA/" + marker_files[marker_ind],header=None).iloc[:, 0].values)
    
# Convert human genes to mouse genes
for pathway in annotations["global"]:
    annotations["global"][pathway] = [g.capitalize() for g in annotations["global"][pathway]]
    
# Only keep HSC and global gene sets
annotations = {key: annotations[key] for key in ["HSC", "global"]}

# Filter gene sets
annotations = check_gene_set_dictionary(adata, annotations, obs_key='cell_type',global_key='global')
for cell_type in ["HSC", "global"]:
    gene_sets = copy.deepcopy(list(annotations[cell_type].keys()))
    for gene_set in gene_sets:
        if len(annotations[cell_type][gene_set]) < 5:
            del annotations[cell_type][gene_set]          
            
# Model fitting
print("Fitting model")
#fit the model (We will run this with only 2 epochs to decrease runtime in this tutorial)
model = spc.est_spectra(adata = adata, gene_set_dictionary = annotations, 
                        use_highly_variable = False, cell_type_key = "cell_type", 
                        use_weights = True, lam = 10, 
                        delta=0.001,kappa = 0.00001, rho = 0.00001, 
                        use_cell_types = True, n_top_vals = 50, 
                        label_factors = True, #whether to label the factors by their overlap coefficient with the input gene sets
                        overlap_threshold = 0.2, #minimum overlap coefficient that has to be surpassed to assign a label to a factor
                        num_epochs=2000 #for demonstration purposes we will only run 2 epochs, we recommend 10,000 epochs
                       )

# Save model to file
print("Saving model to file")
adata.write_h5ad("../../data/mHSCAging10xMultiome/spectra/spectra_adata.h5ad")

# Get markers
factor_markers = adata.uns['SPECTRA_markers']
factor_markers = pd.DataFrame(factor_markers).T
factor_markers.columns = adata.uns['SPECTRA_overlap'].index

# Write markers to file
path = "../../data/mHSCAging10xMultiome/spectra/spectra_markers_lam_10.tsv"
factor_markers.to_csv(path, sep = "\t", index=False)