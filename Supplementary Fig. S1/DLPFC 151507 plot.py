import scanpy as sc
import numpy as np
import pandas as pd
from preprocessing import preprocessing_data
import matplotlib.pyplot as plt
import seaborn as sns


adata = sc.read_visium('WARGA/DLPFC dataset/151507')
adata = preprocessing_data(adata, n_top_genes=2000)
# load ground_truth
ground_truth = pd.read_csv('SpatialDE_clustering/cluster_labels_151507.csv', index_col=0)
true_label = list(ground_truth['ground_truth'])


bayes = pd.read_csv('BayesSpace/BayesSpace_151507_pred_label.csv', index_col=0)
bayes_label = list(bayes['x'])
bayes_label = [x-1 for x in bayes_label]
bayes_label = [str(x) for x in bayes_label]

deepst = sc.read_h5ad('deepST/151507_deepst_results.h5ad')
deepst_label = list(deepst.obs['DeepST_refine_domain'])
deepst_label = [str(x) for x in deepst_label]


sedr = pd.read_csv('SEDR/DLPFC/151507/SEDR/metadata.tsv', index_col=0, sep='\t')
sedr_label = list(sedr['SEDR'])
sedr_label = [str(x) for x in sedr_label]


spaceflow = sc.read_h5ad('spaceflow/151507_spaceflow_results.h5ad')
spaceflow_label = list(spaceflow.obs['spaceflow_pred_label'])
spaceflow_label = [str(x) for x in spaceflow_label]


spaGCN = sc.read_h5ad('spaGCN/DLPFC/151507/SpaGCN/results.h5ad')
spagcn_label = list(spaGCN.obs['refined_pred'])
spagcn_label = [str(x) for x in spagcn_label]


stlearn = pd.read_csv('stlearn/DLPFC/151507/stlearn/metadata.tsv', sep='\t')
stlearn_label = list(stlearn['X_pca_kmeans'])
stlearn_label = [str(x) for x in stlearn_label]


warga = sc.read_h5ad('WARGA/DLPFC dataset/151507/PearlST_results.h5ad')
warga_label = warga.obs['pred_label']
warga_label = [str(x) for x in warga_label]


scanpy = sc.read_h5ad('scanpy/151507_scanpy_results.h5ad')
scanpy_label = list(scanpy.obs['leiden'])
scanpy_label = [str(x) for x in scanpy_label]


seurat = pd.read_csv('Seurat/DLPFC/151507/Seurat/metadata.tsv', sep='\t')
seurat_label = list(seurat['seurat_clusters'])
seurat_label = [str(x) for x in seurat_label]


stagate = sc.read_h5ad('STAGATE/151507_STAGATE_tensorflow_results.h5ad')
stagate_label = list(stagate.obs['mclust'])
stagate_label = [str(x) for x in stagate_label]


specmix = pd.read_csv('SpecMix/151507/pred_label_151507.txt')
specmix_label = list(specmix['label SpiceMixPlus'])
specmix_label = [str(x) for x in specmix_label]



from sklearn.metrics import adjusted_rand_score
adata.obs['Manual annotation'] = true_label
adata.obs['PearlST'] = warga_label
adata.obs['bayes'] = bayes_label
adata.obs['deepst'] = deepst_label
adata.obs['sedr'] = sedr_label
adata.obs['spaceflow'] = spaceflow_label
adata.obs['stlearn'] = stlearn_label
adata.obs['spagcn'] = spagcn_label
adata.obs['scanpy'] = scanpy_label
adata.obs['seurat'] = seurat_label
adata.obs['stagate'] = stagate_label
adata.obs['specmix'] = specmix_label
adata_obs = adata.obs.dropna()

print("PearlST:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['PearlST']))
print("bayes:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['bayes']))
print("deepst:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['deepst']))
print("sedr:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['sedr']))
print("spaceflow:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['spaceflow']))
print("stlearn:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['stlearn']))
print("spagcn:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['spagcn']))
print("scanpy:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['scanpy']))
print("seurat:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['seurat']))
print("stagate:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['stagate']))
print("specmix:  ARI:", adjusted_rand_score(adata_obs['Manual annotation'], adata_obs['specmix']))


plot_true_label_color = ['#A1A9D0', '#F0988C', '#B883D4', '#CFEAF1', '#C4A5DE', '#F6CAE5', '#96CCCB', '#9E9E9E']
adata.uns['Manual annotation_colors'] = plot_true_label_color
plot_color = ['#A1A9D0', '#F0988C', '#B883D4', '#CFEAF1', '#C4A5DE', '#F6CAE5', '#96CCCB']
adata.uns['PearlST_colors'] = plot_color
adata.uns['bayes_colors'] = plot_color
adata.uns['deepst_colors'] = plot_color
adata.uns['sedr_colors'] = plot_color
adata.uns['spaceflow_colors'] = plot_color
adata.uns['stlearn_colors'] = plot_color
adata.uns['spagcn_colors'] = plot_color
adata.uns['scanpy_colors'] = plot_color
adata.uns['seurat_colors'] = plot_color
adata.uns['stagate_colors'] = plot_color
adata.uns['specmix_colors'] = plot_color


fig, axs = plt.subplots(2, 6, figsize=(13,4),constrained_layout=True)
sc.pl.spatial(adata, color=['Manual annotation'], ax=axs[0,0], show=False, spot_size=150, img_key=None)
axs[0,0].spines['right'].set_visible(False) # 去掉边框
axs[0,0].spines['top'].set_visible(False)   # 去掉边框
axs[0,0].spines['left'].set_visible(False) # 去掉边框
axs[0,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('Manual annotation', fontsize=10)


sc.pl.spatial(adata, img_key=None, color=['sedr'], show=False, spot_size=150, ax=axs[0,1],
              legend_fontsize=0, legend_loc='on data')
axs[0,1].spines['right'].set_visible(False) # 去掉边框
axs[0,1].spines['top'].set_visible(False)   # 去掉边框
axs[0,1].spines['left'].set_visible(False) # 去掉边框
axs[0,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('SEDR  ARI:0.429', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['bayes'], show=False, spot_size=150, ax=axs[0,2],
              legend_fontsize=0, legend_loc='on data')
axs[0,2].spines['right'].set_visible(False) # 去掉边框
axs[0,2].spines['top'].set_visible(False)   # 去掉边框
axs[0,2].spines['left'].set_visible(False) # 去掉边框
axs[0,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,2].get_yaxis().set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].set_title('BayesSpace  ARI:0.468', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['spaceflow'], show=False, spot_size=150, ax=axs[0,3],
              legend_fontsize=0, legend_loc='on data')
axs[0,3].spines['right'].set_visible(False) # 去掉边框
axs[0,3].spines['top'].set_visible(False)   # 去掉边框
axs[0,3].spines['left'].set_visible(False) # 去掉边框
axs[0,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,3].get_yaxis().set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].set_title('SpaceFlow  ARI:0.477', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['deepst'], show=False, spot_size=150, ax=axs[0,4],
              legend_fontsize=0, legend_loc='on data')
axs[0,4].spines['right'].set_visible(False) # 去掉边框
axs[0,4].spines['top'].set_visible(False)   # 去掉边框
axs[0,4].spines['left'].set_visible(False) # 去掉边框
axs[0,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,4].get_yaxis().set_visible(False)
axs[0,4].get_xaxis().set_visible(False)
axs[0,4].set_title('DeepST  ARI:0.610', fontsize=10)


sc.pl.spatial(adata, img_key=None, color=['spagcn'], show=False, spot_size=150, ax=axs[0,5],
              legend_fontsize=0, legend_loc='on data')
axs[0,5].spines['right'].set_visible(False) # 去掉边框
axs[0,5].spines['top'].set_visible(False)   # 去掉边框
axs[0,5].spines['left'].set_visible(False) # 去掉边框
axs[0,5].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,5].get_yaxis().set_visible(False)
axs[0,5].get_xaxis().set_visible(False)
axs[0,5].set_title('SpaGCN  ARI:0.475', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['PearlST'], show=False, spot_size=150, ax=axs[1,0],
              )
axs[1,0].spines['right'].set_visible(False) # 去掉边框
axs[1,0].spines['top'].set_visible(False)   # 去掉边框
axs[1,0].spines['left'].set_visible(False) # 去掉边框
axs[1,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('PearlST  ARI:0.596', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['seurat'], show=False, spot_size=150, ax=axs[1,1],
              legend_fontsize=0, legend_loc='on data')
axs[1,1].spines['right'].set_visible(False) # 去掉边框
axs[1,1].spines['top'].set_visible(False)   # 去掉边框
axs[1,1].spines['left'].set_visible(False) # 去掉边框
axs[1,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('Seurat  ARI:0.326', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['stagate'], show=False, spot_size=150, ax=axs[1,2],
              legend_fontsize=0, legend_loc='on data')
axs[1,2].spines['right'].set_visible(False) # 去掉边框
axs[1,2].spines['top'].set_visible(False)   # 去掉边框
axs[1,2].spines['left'].set_visible(False) # 去掉边框
axs[1,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,2].get_yaxis().set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].set_title('STAGATE  ARI:0.210', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['scanpy'], show=False, spot_size=150, ax=axs[1,3],
              legend_fontsize=0, legend_loc='on data')
axs[1,3].spines['right'].set_visible(False) # 去掉边框
axs[1,3].spines['top'].set_visible(False)   # 去掉边框
axs[1,3].spines['left'].set_visible(False) # 去掉边框
axs[1,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,3].get_yaxis().set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].set_title('Scanpy  ARI:0.318', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['specmix'], show=False, spot_size=150, ax=axs[1,4],
              legend_fontsize=0, legend_loc='on data')
axs[1,4].spines['right'].set_visible(False) # 去掉边框
axs[1,4].spines['top'].set_visible(False)   # 去掉边框
axs[1,4].spines['left'].set_visible(False) # 去掉边框
axs[1,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,4].get_yaxis().set_visible(False)
axs[1,4].get_xaxis().set_visible(False)
axs[1,4].set_title('SpecMix  ARI:0.477', fontsize=10)

sc.pl.spatial(adata, img_key=None, color=['stlearn'], show=False, spot_size=150, ax=axs[1,5],
              legend_fontsize=0, legend_loc='on data')
axs[1,5].spines['right'].set_visible(False) # 去掉边框
axs[1,5].spines['top'].set_visible(False)   # 去掉边框
axs[1,5].spines['left'].set_visible(False) # 去掉边框
axs[1,5].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,5].get_yaxis().set_visible(False)
axs[1,5].get_xaxis().set_visible(False)
axs[1,5].set_title('stLearn  ARI:0.493', fontsize=10)


plt.savefig('151507_benchmarking.pdf', dpi=300)



