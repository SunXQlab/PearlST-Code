import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM, plot_pSM_multi


adata = sc.read_h5ad('WARGA/stereoseq_used_barcode_data/stereoseq_used_barcode_epoch_70_results.h5ad')
kmeans = KMeans(n_clusters=7, random_state=0).fit(adata.obsm["WARGA_embed"])
predict_labels = kmeans.predict(adata.obsm["WARGA_embed"])
kmeans_sil = silhouette_score(adata.obsm["WARGA_embed"], predict_labels)
print("kmenas-sil=", "{:.5f}".format(kmeans_sil))

labels = list(predict_labels)
labels = [str(x) for x in labels]
adata.obs['WARGA_refine_domain'] = labels
sc.pl.embedding(adata, basis='spatial', color='WARGA_refine_domain', size=15)


adata_stagate = sc.read_h5ad('STAGATE/stereoseq/stereoseq_used_barcode_STAGATE_tensorflow_results.h5ad')
adata_spaceflow = sc.read_h5ad('spaceflow/stereoseq/stereoseq_spaceflow_results.h5ad')


labels = list(adata.obs['WARGA_refine_domain'])
for i in range(len(labels)):
    if labels[i] == '3':
        labels[i] = 'RMS'
    elif labels[i] == '2':
        labels[i] = 'GCL'
    elif labels[i] == '4':
        labels[i] = 'IPL'
    elif labels[i] == '0':
        labels[i] = 'MCL'
    elif labels[i] == '1':
        labels[i] = 'EPL'
    elif labels[i] == '6':
        labels[i] = 'GL'
    elif labels[i] == '5':
        labels[i] = 'ONL'
adata.obs['WARGA_refine_domain'] = labels
plot_color = ['#A1A9D0', '#F0988C', '#B883D4', '#CFEAF1', '#C4A5DE', '#F6CAE5', '#96CCCB']
# plot_color = ['#F27970', '#FA7F6F', '#54B354', '#05B9E2', '#8983BF', '#C76DA2', '#A1A9D0']
adata.uns['WARGA_refine_domain_colors'] = plot_color
adata_stagate.uns['mclust_colors'] = plot_color
adata_spaceflow.uns['spacefloe_pred_label_colors'] = plot_color

'''
sc.pl.embedding(adata, basis='spatial', color='WARGA_refine_domain', size=15)
sc.pl.embedding(adata_stagate, basis='spatial', color='mclust', size=15)
sc.pl.embedding(adata_spaceflow, basis='spatial', color='spacefloe_pred_label', size=15)
'''

fig, axs = plt.subplots(1, 3, figsize=(11,2.5),constrained_layout=True)
sc.pl.embedding(adata_stagate, basis="spatial",  color='mclust', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('STAGATE', fontsize=12)

sc.pl.embedding(adata_spaceflow, basis="spatial",  color='spacefloe_pred_label', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaceFlow', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='WARGA_refine_domain', size=12, ax=axs[2])
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('PearlST', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/stereoseq_usedbarcode_1.pdf', dpi=300)


fig, axs = plt.subplots(1, 7, figsize=(15,1.8),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='RMS',  color='WARGA_refine_domain', size=12, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('RMS', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='GCL',  color='WARGA_refine_domain', size=12, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('GCL', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='IPL',  color='WARGA_refine_domain', size=12, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('IPL', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='MCL',  color='WARGA_refine_domain', size=12, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('MCL', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='EPL',  color='WARGA_refine_domain', size=12, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('EPL', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='GL',  color='WARGA_refine_domain', size=12, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('GL', fontsize=12)


sc.pl.embedding(adata, basis="spatial", groups='ONL',  color='WARGA_refine_domain', size=12, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].set_title('ONL', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/stereoseq_usedbarcode_2.pdf', dpi=300)



fig, axs = plt.subplots(1, 7, figsize=(15,1.8),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", color='Mbp', size=15, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Mbp', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Nrgn', size=15, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Nrgn', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Pcp4', size=15, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Pcp4', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Gabra1', size=15, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Gabra1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Doc2g', size=15, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('Doc2g', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Cck', size=15, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('S100a5', fontsize=12)


sc.pl.embedding(adata, basis="spatial", color='Apod', size=15, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].set_title('Apod', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/stereoseq_usedbarcode_3.pdf', dpi=300)


pseudo_Spatiotemporal_Map(adata, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_stagate, emb_name='STAGATE', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_spaceflow, emb_name='spaceflow_emb', n_neighbors=20, resolution=1.0)
plot_pSM_multi(adata=adata_stagate, adata1=adata_spaceflow, adata2=adata, adata3=None, scatter_sz=10., rsz=3.,
             csz=11, wspace=.4, hspace=.5, left=0.008, right=0.98, bottom=0.1, top=0.9)


plot_pSM(adata, scatter_sz=5, rsz=5.,
             csz=7., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9)

plt.savefig('C:/Users/haiyu/Desktop/stereo pSM_legend.pdf', dpi=300)





























