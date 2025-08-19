import scanpy as sc 
import os
import anndata as ad
import pandas as pd
import pickle

sample_names = ["PBMC_10k"]
data_paths = ["/gstore/project/epigen/PBMC/raw_data/"]
barcode_tab = pd.read_csv("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/barcode_tab_PBMC.csv")
adata_list = list(map(lambda x: sc.read_10x_h5(os.path.join(x, "pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")), data_paths))
for i in range(len(sample_names)):
    adata_list[i].obs['sample_id'] = sample_names[i]
    adata_list[i].obs_names = list(map(lambda x: x.replace("-", "."), adata_list[i].obs_names))
    sample_barcodes = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]]["barcode"]
    if False in list(map(lambda x: x in adata_list[i].obs_names, sample_barcodes)):
        raise Exception("All barcodes should be present in the gene expression data")
    adata_list[i] = adata_list[i][sample_barcodes,]
    adata_list[i].obs["cell_type"] = barcode_tab.loc[barcode_tab['sample_id'] == sample_names[i]]["cell_type"].values
    adata_list[i].obs_names = list(map(lambda x: x[0]+"___"+x[1], zip(adata_list[i].obs_names, adata_list[i].obs['sample_id'])))
    adata_list[i].var_names_make_unique()
adatas = dict(zip(sample_names, adata_list))
adata = ad.concat(adatas, label = "dataset")
adata.var_names_make_unique()
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write_h5ad("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/adata.h5ad")
