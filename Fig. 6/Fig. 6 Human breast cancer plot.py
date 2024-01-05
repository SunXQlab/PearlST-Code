"""
plot benchmarking Human_breast_cancer
"""


import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib.pyplot as plt
from utils import pseudo_Spatiotemporal_Map, plot_pSM
from preprocessing import preprocessing_data

adata = sc.read_visium('./WARGA/DLPFC dataset/Human breast cancer')
adata = preprocessing_data(adata, n_top_genes=2000)
true_label = pd.read_csv('./WARGA/DLPFC dataset/Human breast cancer/metadata.tsv', sep='\t', index_col=0)
adata.obs['Manual annotation'] = list(true_label['ground_truth'])

bayes = pd.read_csv('BayesSpace/BayesSpace_human_breast_cancer_pred_label.csv', index_col=0)
bayes_label = list(bayes['x'])
bayes_label = [x-1 for x in bayes_label]
bayes_label = [str(x) for x in bayes_label]

deepst = sc.read_h5ad('deepST/Human breast cancer_deepst_results.h5ad')
deepst_label = list(deepst.obs['DeepST_refine_domain'])
deepst_label = [str(x) for x in deepst_label]

sedr = sc.read_h5ad('SEDR/DLPFC/human breast cancer/human_breast_cancer_SEDR_results.h5ad')
sedr_label = list(sedr.obs['SEDR_label'])
sedr_label = [str(x) for x in sedr_label]

spaceflow = sc.read_h5ad('spaceflow/human_breast_cancer_spaceflow_results.h5ad')
spaceflow_label = list(spaceflow.obs['spaceflow_pred_label'])
spaceflow_label = [str(x) for x in spaceflow_label]

spaGCN = sc.read_h5ad('spaGCN/DLPFC/Human breast cancer/SpaGCN/results.h5ad')
spagcn_label = list(spaGCN.obs['refined_pred'])
spagcn_label = [str(x) for x in spagcn_label]

stlearn = pd.read_csv('stlearn/DLPFC/Human breast cancer/stlearn/metadata.tsv', sep='\t')
stlearn_label = list(stlearn['X_pca_kmeans'])
stlearn_label = [str(x) for x in stlearn_label]

warga = sc.read_h5ad('WARGA/DLPFC dataset/Human breast cancer/PearlST_results.h5ad')
warga_label = warga.obs['pred_label']
warga_label = [str(x) for x in warga_label]

scanpy = sc.read_h5ad('scanpy/Human_breast_cancer_scanpy_results.h5ad')
scanpy_label = list(scanpy.obs['leiden'])
scanpy_label = [str(x) for x in scanpy_label]

seurat = pd.read_csv('Seurat/DLPFC/Human breast cancer/metadata.tsv', sep='\t')
seurat_label = list(seurat['seurat_clusters'])
seurat_label = [str(x) for x in seurat_label]

stagate = sc.read_h5ad('STAGATE/human_breast_cancer_STAGATE_tensorflow_results.h5ad')
stagate_label = list(stagate.obs['mclust'])
stagate_label = [str(x) for x in stagate_label]



from sklearn.metrics import adjusted_rand_score
adata.obs['PearlST'] = warga_label
# adata.obs['bayes'] = bayes_label
adata.obs['deepst'] = deepst_label
adata.obs['sedr'] = sedr_label
adata.obs['spaceflow'] = spaceflow_label
adata.obs['stlearn'] = stlearn_label
adata.obs['spagcn'] = spagcn_label
adata.obs['scanpy'] = scanpy_label
adata.obs['seurat'] = seurat_label
adata.obs['stagate'] = stagate_label
# adata.obs['specmix'] = specmix_label
adata_obs = adata.obs.dropna()

print("PearlST:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['PearlST']))
# print("bayes:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['bayes']))
print("deepst:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['deepst']))
print("sedr:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['sedr']))
print("spaceflow:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['spaceflow']))
print("stlearn:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['stlearn']))
print("spagcn:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['spagcn']))
print("scanpy:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['scanpy']))
print("seurat:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['seurat']))
print("stagate:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['stagate']))
# print("specmix:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['specmix']))



fig, axs = plt.subplots(1, 2, figsize=(13,4),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='Manual annotation', size=50, ax=axs[0])
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Manual annotation', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='PearlST', size=50, ax=axs[1])
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('PearlST  ARI:0.653', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/human breast cancer_1.pdf', dpi=300)



fig, axs = plt.subplots(1, 5, figsize=(13,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", color='deepst', size=30, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('DeepST  ARI:0.559', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='sedr', size=30, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SEDR  ARI:0.494', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='spaceflow', size=30, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('SpaceFlow  ARI:0.491', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='stlearn', size=30, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].invert_yaxis()
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('stLearn  ARI:0.541', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='spagcn', size=30, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].invert_yaxis()
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('SpaGCN  ARI:0.456', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/human breast cancer_2.pdf', dpi=300)



fig, axs = plt.subplots(1, 3, figsize=(7.8,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", color='scanpy', size=30, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].invert_yaxis()
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Scanpy  ARI:0.496', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='seurat', size=30, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].invert_yaxis()
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Seurat  ARI:0.474', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='stagate', size=30, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].invert_yaxis()
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('STAGATE  ARI:0.343', fontsize=12)


plt.savefig('C:/Users/haiyu/Desktop/human breast cancer_3.pdf', dpi=300)


