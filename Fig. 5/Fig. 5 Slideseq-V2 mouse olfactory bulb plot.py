"""
plot STAGATE slideseq-V2
"""


import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt


adata = sc.read_h5ad('WARGA/slideseq-V2 STAGATE/Slideseq-V2_Puck_epoch300_results.h5ad')

kmeans = KMeans(n_clusters=11, random_state=0).fit(adata.obsm["WARGA_embed"])
predict_labels = kmeans.predict(adata.obsm["WARGA_embed"])
adata.obs['pred_label'] = [str(x) for x in predict_labels]
kmeans_sil = silhouette_score(adata.obsm["WARGA_embed"], adata.obs['pred_label'])
print("kmenas-sil=", "{:.5f}".format(kmeans_sil))



sc.pl.embedding(adata, basis='spatial', color='pred_label', size=10)
ax = plt.gca()
ax.invert_yaxis()


labels = list(adata.obs['pred_label'])
for i in range(len(labels)):
    if labels[i] == '0':
        labels[i] = 'RMS'
    elif labels[i] == '1':
        labels[i] = 'GL_2'
    elif labels[i] == '2':
        labels[i] = 'EPL'
    elif labels[i] == '3':
        labels[i] = 'AOBgr'
    elif labels[i] == '4':
        labels[i] = 'MCL'
    elif labels[i] == '5':
        labels[i] = 'GCL_2'
    elif labels[i] == '6':
        labels[i] = 'ONL'
    elif labels[i] == '7':
        labels[i] = 'GCL_1'
    elif labels[i] == '8':
        labels[i] = 'AOB'
    elif labels[i] == '9':
        labels[i] = 'GL_1'
    elif labels[i] == '10':
        labels[i] = 'GCL_1'
adata.obs['pred_label'] = labels


sc.pl.embedding(adata, basis='spatial', color='pred_label', size=10)
ax = plt.gca()
ax.invert_yaxis()
ax.spines['right'].set_visible(False) # 去掉边框
ax.spines['top'].set_visible(False)   # 去掉边框
ax.spines['left'].set_visible(False) # 去掉边框
ax.spines['bottom'].set_visible(False)   # 去掉边框
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)
ax.set_title('PearlST annotation', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 Puck PearlST annotation.pdf', dpi=300)



fig, axs = plt.subplots(1, 7, figsize=(13,2),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='MCL',  color='pred_label', size=5, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('MCL', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='GCL_2',  color='pred_label', size=5, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('GCL_2', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='RMS',  color='pred_label', size=5, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('RMS', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='AOBgr',  color='pred_label', size=5, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('AOBgr', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='AOB',  color='pred_label', size=5, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('AOB', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='EPL',  color='pred_label', size=5, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].invert_yaxis()
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('EPL', fontsize=12)


sc.pl.embedding(adata, basis="spatial", groups='GCL_1',  color='pred_label', size=8, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].invert_yaxis()
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].set_title('GCL_1', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 Puck groups.pdf', dpi=300)



fig, axs = plt.subplots(1, 7, figsize=(13,1.8),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", color='Gabra1', size=5, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Gabra1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Slc6a11', size=5, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Slc6a11', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='Mbp', size=5, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Mbp', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='Atp2b4', size=5, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Atp2b4', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='Gap43', size=5, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('Gap43', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='Pcp4', size=5, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].invert_yaxis()
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].set_title('Pcp4', fontsize=12)


sc.pl.embedding(adata, basis="spatial",  color='Nrgn', size=5, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].invert_yaxis()
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].set_title('Nrgn', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/Slideseq-V2 Puck groups_markers.pdf', dpi=300)
