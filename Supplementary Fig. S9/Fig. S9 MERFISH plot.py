import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM


adata = sc.read_h5ad('WARGA/MERFISH/MERFISH_Allen2022Molecular_results.h5ad')

fig, axs = plt.subplots(1, 3, figsize=(13,4),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='cell_type_annot', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Cell type annotation', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='ori_pred_label', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Raw data (ARI: 0.7222)', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='augment_pred_label', size=12, ax=axs[2])
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Augment data (ARI: 0.8031)', fontsize=14)

plt.savefig('MERFISH_Allen2022Molecular.pdf', dpi=300)



adata = sc.read_h5ad('WARGA/MERFISH/MERFISH_moffitt2018molecular_results.h5ad')

fig, axs = plt.subplots(1, 3, figsize=(13,4),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='Cell_class', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Cell type annotation', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='ori_pred_label', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Raw data (ARI: 0.3717)', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='augment_pred_label', size=12, ax=axs[2])
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Augment data (ARI: 0.4201)', fontsize=14)

plt.savefig('MERFISH_moffitt2018molecular.pdf', dpi=300)


