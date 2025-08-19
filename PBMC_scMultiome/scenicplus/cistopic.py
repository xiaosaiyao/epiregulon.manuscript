import scanpy as sc 
import os
import anndata as ad
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
import dill
import warnings
import pandas as pd


import requests
import numpy as np
import pybiomart as pbm
import pickle
import ray

n_cores = int(sys.argv[1])

sample_names = ["PBMC_10k"]
adata = sc.read_h5ad("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/adata.h5ad")
paths_to_fragments = ["/gstore/project/epigen/PBMC/raw_data/"]
paths_to_peak_matrix = ["/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/peak_matrix.tsv"]

scRNA_bc = adata.obs_names
cell_data = adata.obs
cell_data['cell_type'] = cell_data['cell_type'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.
fragments_dict = dict(zip(sample_names, paths_to_fragments))
matrices_dict = dict(zip(sample_names, paths_to_peak_matrix))
cistopic_obj_list = [create_cistopic_object_from_matrix_file(fragment_matrix_file = matrices_dict[key],
                                                  path_to_fragments=fragments_dict[key],
                                                  project = key) for key in fragments_dict.keys()]

cistopic_obj = merge(cistopic_obj_list)
cistopic_obj.add_cell_data(cell_data)
f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/cistopic_obj_PBMC", "wb")
pickle.dump(cistopic_obj, f)
f.close()
models=run_cgs_models(cistopic_obj,
            n_topics=[25],
            n_cpu=n_cores,
            n_iter=500,
            random_state=555,
            alpha=50,
            alpha_by_topic=True,
            eta=0.1,
            eta_by_topic=False,
            save_path=None,
            _temp_dir = '/gstore/project/epigen/temp/scenic/')


f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/models.pkl", "wb")
pickle.dump(models, f)
f.close()

# pycistopic tss get_tss  --output /gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/qc/tss.bed --name "hsapiens_gene_ensembl" --to-chrom-source ucsc --ucsc hg38


f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/models.pkl", "rb")
models = pickle.load(f)
f.close()


f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/cistopic_obj_PBMC", "rb")
cistopic_obj = pickle.load(f)
f.close()


# results of the model evalution available in /gstore/project/epigen/PBMC/analysis/scenicplus/Model_evaluation.ipynb

model = evaluate_models(models,
                   select_model=25,
                   return_model=True,
                   metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                   plot_metrics=False)

cistopic_obj.add_LDA_model(model)

f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/cistopic_obj_PBMC", "wb")
pickle.dump(cistopic_obj, f)
f.close()
