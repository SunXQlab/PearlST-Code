import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM, plot_pSM_multi


adata = sc.read_h5ad('WARGA/STARmap/STARmap_results.h5ad')
ari_adata = adjusted_rand_score(adata.obs['label'], adata.obs['WARGA_refine_domain'])
adata_spaceflow = sc.read_h5ad('spaceflow/STARmap/STAR_spaceflow_results.h5ad')
ari_spaceflow = adjusted_rand_score(adata_spaceflow.obs['label'], adata_spaceflow.obs['spacefloe_pred_label'])
adata_stagate = sc.read_h5ad('STAGATE/STARmap/STARmap_STAGATE_tensorflow_results.h5ad')
ari_stagate = adjusted_rand_score(adata_stagate.obs['label'], adata_stagate.obs['mclust'])
adata_sedr = sc.read_h5ad('SEDR/STARmap/SEDR_results.h5ad')
ari_sedr = adjusted_rand_score(adata_sedr.obs['label'], adata_sedr.obs['SEDR_label'])
adata_scanpy = sc.read_h5ad('scanpy/STARmap_scanpy_results.h5ad')
ari_scanpy = adjusted_rand_score(adata_scanpy.obs['label'], adata_scanpy.obs['leiden'])


fig, axs = plt.subplots(1, 2, figsize=(13,3.5),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='label', size=50, ax=axs[0])
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Layer structure', fontsize=12)

sc.pl.embedding(adata_spaceflow, basis="spatial",  color='cell_type', size=50, ax=axs[1])
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Cell type', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/STARmap_layers.pdf', dpi=300)



fig, axs = plt.subplots(1, 5, figsize=(14,3),constrained_layout=True)
sc.pl.embedding(adata_spaceflow, basis="spatial", color='spacefloe_pred_label', size=30, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('SpaceFlow (ARI: 0.6051)', fontsize=12)

sc.pl.embedding(adata_stagate, basis="spatial", color='mclust', size=30, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('STAGATE (ARI: 0.2561)', fontsize=12)

sc.pl.embedding(adata_sedr, basis="spatial", color='SEDR_label', size=30, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('SEDR (ARI: 0.3510)', fontsize=12)

sc.pl.embedding(adata_scanpy, basis="spatial", color='leiden', size=30, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Scanpy (ARI: 0.1653)', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='WARGA_refine_domain', size=30, ax=axs[4])
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('PearlST (ARI: 0.6766)', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/STARmap_layer_pred.pdf', dpi=300)

pseudo_Spatiotemporal_Map(adata_spaceflow, emb_name='spaceflow_emb', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_stagate, emb_name='STAGATE', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_sedr, emb_name='SEDR_emb', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_scanpy, emb_name='X_pca', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata, emb_name='WARGA_embed', n_neighbors=10, resolution=1.0)


plot_pSM_multi(adata=adata_spaceflow, adata1=adata_stagate, adata2=adata_sedr, adata3=adata_scanpy, adata4=adata,
               scatter_sz=15, rsz=3., csz=13, wspace=.4, hspace=.5, left=0.008, right=0.98, bottom=0.1, top=0.9
               )




plot_pSM(adata_spaceflow, scatter_sz=20., rsz=5.,
             csz=6.5, wspace=.4, hspace=.5, left=0.088, right=1.03, bottom=0.1, top=0.9)

plt.savefig('C:/Users/haiyu/Desktop/STARmap_pseudo_map_legend.pdf', dpi=300)


