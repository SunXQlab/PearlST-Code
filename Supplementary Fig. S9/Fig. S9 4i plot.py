import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM


adata = sc.read_h5ad('WARGA/4i/4i_119_results.h5ad')
print("raw data:  ARI:", adjusted_rand_score(adata.obs['cluster'], adata.obs['raw_data_label']))
print("raw data:  ARI:", adjusted_rand_score(adata.obs['cluster'], adata.obs['augment_pred_label']))




fig, axs = plt.subplots(1, 3, figsize=(13,3.5),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='cluster', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Ground truth', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='raw_data_label', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Raw data (ARI: 0.558)', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='augment_pred_label', size=12, ax=axs[2])
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Augment data (ARI: 0.596)', fontsize=14)

plt.savefig('4i_119.pdf', dpi=300)



adata = sc.read_h5ad('WARGA/4i/4i_122_results.h5ad')
print("raw data:  ARI:", adjusted_rand_score(adata.obs['cluster'], adata.obs['raw_data_label']))
print("raw data:  ARI:", adjusted_rand_score(adata.obs['cluster'], adata.obs['augment_pred_label']))




fig, axs = plt.subplots(1, 3, figsize=(13,3.5),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='cluster', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Ground truth', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='raw_data_label', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Raw data (ARI: 0.430)', fontsize=14)

sc.pl.embedding(adata, basis="spatial",  color='augment_pred_label', size=12, ax=axs[2])
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Augment data (ARI: 0.474)', fontsize=14)

plt.savefig('4i_122.pdf', dpi=300)

