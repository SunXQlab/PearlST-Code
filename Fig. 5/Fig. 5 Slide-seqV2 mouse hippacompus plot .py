"""
plot GraphST slideseq-V2
"""


import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM


pearst = sc.read_h5ad('WARGA/GraphST/slideseq-V2_GraphST_90_results.h5ad')
label = list(pearst.obs['pred_label'])
label = [str(x) for x in label]
pearst.obs['pred_label'] = label
sil = silhouette_score(pearst.obsm['WARGA_embed'], pearst.obs['pred_label'])

spaceflow = sc.read_h5ad('spaceflow/slideseq-V2/slideseq-V2_GraphST_spaceflow_results.h5ad')
stagate = sc.read_h5ad('STAGATE/slideseq-V2_GraphST_STAGATE_tensorflow_results.h5ad')

'''
plot_color = ['#A1A9D0', '#F0988C', '#B883D4', '#CFEAF1', '#C4A5DE', '#F6CAE5', '#96CCCB', '#82B0D2', '#FFBE7A', '#E7DAD2']
pearst.uns['pred_label_colors'] = plot_color
spaceflow.uns['spacefloe_pred_label_colors'] = plot_color
stagate.uns['mclust_colors'] = plot_color
'''


fig, axs = plt.subplots(1, 3, figsize=(10,2.8),constrained_layout=True)
sc.pl.embedding(stagate, basis="spatial",  color='mclust', size=8, ax=axs[0])
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('STAGATE', fontsize=12)

sc.pl.embedding(spaceflow, basis="spatial",  color='spacefloe_pred_label', size=8, ax=axs[1])
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaceFlow', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", color='pred_label', size=8, ax=axs[2])
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('PearlST', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 GraphST_benchmarking.pdf', dpi=300)



fig, axs = plt.subplots(1, 5, figsize=(13,2.8),constrained_layout=True)
sc.pl.embedding(pearst, basis="spatial", groups='4',  color='pred_label', size=8, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('CA1', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", groups='1',  color='pred_label', size=8, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('CA3', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", groups='7',  color='pred_label', size=8, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('DG', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", groups='8',  color='pred_label', size=8, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('V3', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", groups='9',  color='pred_label', size=8, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('LH', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 GraphST_groups.pdf', dpi=300)




fig, axs = plt.subplots(1, 5, figsize=(13,2.5),constrained_layout=True)
sc.pl.embedding(pearst, basis="spatial", color='Wfs1', size=8, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Wfs1', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", color='Cpne4', size=8, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Cpne4', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", color='C1ql2', size=8, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('C1ql2', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", color='Enpp2', size=8, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Enpp2', fontsize=12)

sc.pl.embedding(pearst, basis="spatial", color='Nwd2', size=8, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('Nwd2', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 GraphST_groups_markers.pdf', dpi=300)







pseudo_Spatiotemporal_Map(spaceflow, pSM_values_save_filepath="./pSM_values.tsv", n_neighbors=20, resolution=1.0)
plot_pSM(spaceflow, pSM_figure_save_filepath="./pseudo-Spatiotemporal-Map.pdf", colormap='roma', scatter_sz=1., rsz=5.,
             csz=7., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

plt.savefig('Slideseq-V2 GraphST_SpaceFlow_pseudo_map.pdf', dpi=300)


