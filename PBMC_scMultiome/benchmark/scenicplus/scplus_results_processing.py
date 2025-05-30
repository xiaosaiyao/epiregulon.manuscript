import mudata
import pandas as pd
scplus_mdata = mudata.read("/gstore/project/epigen/benchmark/PBMC/scenicplus/scplus_pipeline/Snakemake/scplusmdata.h5mu")

colnames = scplus_mdata['direct_gene_based_AUC'].var_names.str.replace("+","pos").str.replace("-","neg")
pd.DataFrame(scplus_mdata['direct_gene_based_AUC'].X, index = scplus_mdata['direct_gene_based_AUC'].obs_names, columns = colnames).to_csv('/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/gene_based_AUC.csv')
colnames = scplus_mdata['direct_region_based_AUC'].var_names.str.replace("+","pos").str.replace("-","neg")
pd.DataFrame(scplus_mdata['direct_region_based_AUC'].X, index = scplus_mdata['direct_region_based_AUC'].obs_names, columns = colnames).to_csv('/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/region_based_AUC.csv')
scplus_mdata.uns['direct_e_regulon_metadata'].to_csv('/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/scenic_GRN.csv')
scplus_mdata.uns['extended_e_regulon_metadata'].to_csv('/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/scenic_GRN_extended.csv')
