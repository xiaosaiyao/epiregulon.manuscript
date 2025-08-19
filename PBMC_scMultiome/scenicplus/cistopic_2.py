from pycistarget.utils import region_names_to_coordinates
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.clust_vis import *
from pycisTopic.topic_binarization import *
from pycisTopic.diff_features import *
import pyranges as pr
import os
import pandas as pd

n_cores = int(sys.argv[1])

f = open("/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/cistopic_obj_PBMC", "rb")
cistopic_obj = pickle.load(f)
f.close()


run_umap(cistopic_obj, target  = 'cell', scale=True)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=False, num_columns=5
)


region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3000,
    plot=False, num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False, min_disp=0.05,
min_mean = 0.0125, max_mean = 3, max_disp = np.inf, n_bins =20, n_top_features = None)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='cell_type', var_features=variable_regions,
contrasts = None, adjpval_thr = 0.05, log2fc_thr=np.log2(1.5), n_cpu=n_cores,  _temp_dir='/gstore/project/epigen/temp/',
    split_pattern = '-')
out_dir = "/gstore/project/epigen/benchmark/PBMC/OUTPUT/scenicplus/"

os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
    
for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity

chromsizes = pd.read_table(os.path.join(out_dir, "qc", "hg38.chrom_sizes_and_alias.tsv"))
chromsizes.rename({"# ucsc": "Chromosome", "length": "End"}, axis = 1, inplace = True)
chromsizes["Start"] = 0
chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])

pr_annotation = pd.read_table(
        os.path.join(out_dir, "qc", "tss.bed")
    ).rename(
        {"Name": "Gene", "# Chromosome": "Chromosome"}, axis = 1)
pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
pr_annotation = pr.PyRanges(pr_annotation)


gene_act, weigths = get_gene_activity(
    imputed_acc_obj,
    pr_annotation,
    chromsizes,
    use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
    upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                             # these bp will be taken (1kbp here)
    downstream=[1000,100000], # Search space downstream
    distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
    decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
    extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                          #this weight)
    extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
    gene_size_weight=False, # Whether to add a weights based on the length of the gene
    gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                          #in the genome
    remove_promoters=False, # Whether to remove promoters when computing gene activity scores
    average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                          #activity score
    scale_factor=1, # Value to multiply for the final gene activity matrix
    extend_tss=[10,10], # Space to consider a promoter
    gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
    return_weights= True, # Whether to return the final weights
    project='Gene_activity') # Project name for the gene activity object
    
DAG_markers_dict= find_diff_features(
    cistopic_obj,
    gene_act,
    variable='cell_type',
    var_features=None,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=n_cores,
    _temp_dir='/gstore/scratch/u/wlodarct/ray_spill',
    split_pattern = '-')

from pycisTopic.loom import export_region_accessibility_to_loom, export_gene_activity_to_loom
cluster_markers = {'cell_type': markers_dict}
os.makedirs(os.path.join(out_dir, "loom"), exist_ok=True)

export_region_accessibility_to_loom(
    accessibility_matrix = imputed_acc_obj,
    cistopic_obj = cistopic_obj,
    binarized_topic_region = region_bin_topics_otsu,
    binarized_cell_topic = binarized_cell_topic,
    selected_cells = cistopic_obj.projections['cell']['UMAP'].index.tolist(),
    out_fname = os.path.join(out_dir, "loom", "10x_multiome_blood_pycisTopic_region_accessibility.loom"),
    cluster_annotation = ['cell_type'],
    cluster_markers = cluster_markers,
    tree_structure = ('10x_multiome_blood', 'pycisTopic', 'noDBL_all'),
    title = 'Region accessibility',
    nomenclature = "hg38",
    split_pattern = '-'
)

export_gene_activity_to_loom(
    gene_activity_matrix = gene_act,
    cistopic_obj = cistopic_obj,
    out_fname = os.path.join(out_dir, "loom", "10x_multiome_blood_pycisTopic_gene_activity.loom"),
    cluster_annotation = ['cell_type'],
    cluster_markers = cluster_markers,
    tree_structure = ('10x_multiome_blood', 'pycisTopic', 'ATAC'),
    title = 'Gene activity',
    nomenclature = "hg38",
    split_pattern = '-'
)
