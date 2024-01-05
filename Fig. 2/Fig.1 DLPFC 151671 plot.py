import scanpy as sc
import numpy as np
import pandas as pd
from preprocessing import preprocessing_data
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import adjusted_rand_score


data = {
    "BayesSpace": [0.4680, 0.4374, 0.3944, 0.4341, 0.4711, 0.4297, 0.7335, 0.4391, 0.5494, 0.2957, 0.5289, 0.3674],
    "DeepST": [0.6107, 0.3600, 0.3638, 0.5202, 0.3883, 0.3284, 0.4945, 0.5517, 0.4766, 0.3963, 0.5164, 0.4458],
    "SpaGCN": [0.4755, 0.3669, 0.3769, 0.4433, 0.2240, 0.3725, 0.5458, 0.5606, 0.4608, 0.3234, 0.3004, 0.3257],
    "stLearn": [0.4927, 0.3149, 0.4143, 0.4442, 0.3257, 0.2282, 0.3889, 0.3474, 0.3052, 0.3859, 0.3844, 0.3998],
    "SEDR": [0.4290, 0.3839, 0.3829, 0.3025, 0.3642, 0.4230, 0.4457, 0.4916, 0.5626, 0.4307, 0.4551, 0.4840],
    "SpaceFlow": [0.4770, 0.3938, 0.3938, 0.3496, 0.3140, 0.1585, 0.3268, 0.4154, 0.4276, 0.3160, 0.3800, 0.3142],
    "Seurat": [0.3262, 0.3790, 0.2235, 0.3319, 0.3231, 0.3174, 0.2733, 0.1906, 0.1694, 0.3042, 0.2920, 0.3110],
    "Scanpy": [0.3189, 0.2460, 0.2800, 0.2295, 0.1379, 0.1446, 0.1992, 0.1433, 0.1694, 0.3192, 0.2876, 0.2182],
    "STAGATE": [0.2102, 0.2129, 0.2426, 0.2436, 0.2716, 0.2465, 0.3595, 0.3468, 0.3615, 0.3826, 0.3872, 0.3152],
    "SpecMix": [0.4775, 0.4176, 0.4069, 0.3165, 0.1444, 0.1209, 0.4578, 0.5660, 0.5678, 0.4023, 0.3997, 0.4348],
    "PearlST": [0.5966, 0.5451, 0.5638, 0.5296, 0.7405, 0.7807, 0.8243, 0.7812, 0.6228, 0.5695, 0.5486, 0.5169],
}

df = pd.DataFrame(data)
fig, ax = plt.subplots(figsize=(13,5))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
sns.boxplot(data=df, palette='viridis')
sns.stripplot(data=df, color="orange", jitter=0.2, size=4)
plt.xticks(["BayesSpace", "DeepST", "SpaGCN", "stLearn", "SEDR", "SpaceFlow", "Seurat", "Scanpy", "STAGATE", "SpecMix", "PearlST"],
           rotation=10, fontsize=10)
plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], fontsize=10)
plt.xlabel('Methods', fontsize=12)
plt.ylabel('ARI score', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/true_label_box.pdf', dpi=300)





adata = sc.read_visium('151671')
adata = preprocessing_data(adata, n_top_genes=2000)

# load ground_truth
ground_truth = pd.read_csv('SpatialDE_clustering/cluster_labels_151671.csv', index_col=0)
true_label = list(ground_truth['ground_truth'])
adata.obs['Manual annotation'] = true_label
plot_true_label_color = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC', '#E7DAD2']
adata.uns['Manual annotation_colors'] = plot_true_label_color



bayes = pd.read_csv('BayesSpace/BayesSpace_151671_pred_label.csv', index_col=0)
bayes_label = list(bayes['x'])
bayes_label = [x-1 for x in bayes_label]
bayes_label = [str(x) for x in bayes_label]

for i in range(len(bayes_label)):
    if bayes_label[i] == '4':
        bayes_label[i] = 'WM'
    elif bayes_label[i] == '1':
        bayes_label[i] = 'Layer_6'
    elif bayes_label[i] == '2':
        bayes_label[i] = 'Layer_5'
    elif bayes_label[i] == '0':
        bayes_label[i] = 'Layer_4'
    elif bayes_label[i] == '3':
        bayes_label[i] = 'Layer_3'

deepst_adata = sc.read_h5ad('deepST/151671_deepst.h5ad')
deepst_label = list(deepst_adata.obs['DeepST_refine_domain'])
deepst_label = [str(x) for x in deepst_label]

for i in range(len(deepst_label)):
    if deepst_label[i] == '1':
        deepst_label[i] = 'WM'
    elif deepst_label[i] == '4':
        deepst_label[i] = 'Layer_6'
    elif deepst_label[i] == '2':
        deepst_label[i] = 'Layer_5'
    elif deepst_label[i] == '3':
        deepst_label[i] = 'Layer_4'
    elif deepst_label[i] == '0':
        deepst_label[i] = 'Layer_3'


sedr = pd.read_csv('SEDR/151671/metadata.tsv', index_col=0, sep='\t')
sedr_label = list(sedr['SEDR'])
sedr_label = [str(x) for x in sedr_label]
for i in range(len(sedr_label)):
    if sedr_label[i] == '3':
        sedr_label[i] = 'WM'
    elif sedr_label[i] == '4':
        sedr_label[i] = 'Layer_6'
    elif sedr_label[i] == '2':
        sedr_label[i] = 'Layer_5'
    elif sedr_label[i] == '1':
        sedr_label[i] = 'Layer_4'
    elif sedr_label[i] == '0':
        sedr_label[i] = 'Layer_3'


spaceflow = sc.read_h5ad('spaceflow/151671_spaceflow.h5ad')
spaceflow_label = list(spaceflow.obs['leiden'])
spaceflow_label = [str(x) for x in spaceflow_label]

for i in range(len(sedr_label)):
    if spaceflow_label[i] == '0':
        spaceflow_label[i] = 'WM'
    elif spaceflow_label[i] == '2':
        spaceflow_label[i] = 'Layer_6'
    elif spaceflow_label[i] == '3':
        spaceflow_label[i] = 'Layer_5'
    elif spaceflow_label[i] == '4':
        spaceflow_label[i] = 'Layer_4'
    elif spaceflow_label[i] == '1':
        spaceflow_label[i] = 'Layer_3'

spaGCN = sc.read_h5ad('spaGCN/151671_results.h5ad')
spagcn_label = list(spaGCN.obs['refined_pred'])
spagcn_label = [str(x) for x in spagcn_label]

for i in range(len(sedr_label)):
    if spagcn_label[i] == '4':
        spagcn_label[i] = 'WM'
    elif spagcn_label[i] == '2':
        spagcn_label[i] = 'Layer_6'
    elif spagcn_label[i] == '1':
        spagcn_label[i] = 'Layer_5'
    elif spagcn_label[i] == '3':
        spagcn_label[i] = 'Layer_4'
    elif spagcn_label[i] == '0':
        spagcn_label[i] = 'Layer_3'


stlearn = pd.read_csv('stlearn/151671/metadata.tsv', sep='\t')
stlearn_label = list(stlearn['X_pca_kmeans'])
stlearn_label = [str(x) for x in stlearn_label]

for i in range(len(sedr_label)):
    if stlearn_label[i] == '1':
        stlearn_label[i] = 'WM'
    elif stlearn_label[i] == '2':
        stlearn_label[i] = 'Layer_6'
    elif stlearn_label[i] == '4':
        stlearn_label[i] = 'Layer_5'
    elif stlearn_label[i] == '3':
        stlearn_label[i] = 'Layer_4'
    elif stlearn_label[i] == '0':
        stlearn_label[i] = 'Layer_3'


warga = sc.read_h5ad('WARGA/DLPFC dataset/151671/PearlST_results.h5ad')
warga_label = warga.obs['pred_label']
warga_label = [str(x) for x in warga_label]

for i in range(len(warga_label)):
    if warga_label[i] == '4':
        warga_label[i] = 'WM'
    elif warga_label[i] == '2':
        warga_label[i] = 'Layer_6'
    elif warga_label[i] == '0':
        warga_label[i] = 'Layer_5'
    elif warga_label[i] == '1':
        warga_label[i] = 'Layer_4'
    elif warga_label[i] == '3':
        warga_label[i] = 'Layer_3'


scanpy = sc.read_h5ad('scanpy/151671_scanpy_results.h5ad')
scanpy_label = list(scanpy.obs['leiden'])
scanpy_label = [str(x) for x in scanpy_label]
for i in range(len(scanpy_label)):
    if scanpy_label[i] == '4':
        scanpy_label[i] = 'WM'
    elif scanpy_label[i] == '3':
        scanpy_label[i] = 'Layer_6'
    elif scanpy_label[i] == '2':
        scanpy_label[i] = 'Layer_5'
    elif scanpy_label[i] == '0':
        scanpy_label[i] = 'Layer_4'
    elif scanpy_label[i] == '1':
        scanpy_label[i] = 'Layer_3'

seurat = pd.read_csv('Seurat/DLPFC/151671/Seurat/metadata.tsv', sep='\t')
seurat_label = list(seurat['seurat_clusters'])
seurat_label = [str(x) for x in seurat_label]
for i in range(len(seurat_label)):
    if seurat_label[i] == '4':
        seurat_label[i] = 'WM'
    elif seurat_label[i] == '3':
        seurat_label[i] = 'Layer_6'
    elif seurat_label[i] == '2':
        seurat_label[i] = 'Layer_5'
    elif seurat_label[i] == '0':
        seurat_label[i] = 'Layer_4'
    elif seurat_label[i] == '1':
        seurat_label[i] = 'Layer_3'


stagate = sc.read_h5ad('STAGATE/151671_STAGATE_tensorflow_results.h5ad')
stagate_label = list(stagate.obs['mclust'])
stagate_label = [str(x) for x in stagate_label]
for i in range(len(stagate_label)):
    if stagate_label[i] == '5':
        stagate_label[i] = 'WM'
    elif stagate_label[i] == '2':
        stagate_label[i] = 'Layer_6'
    elif stagate_label[i] == '4':
        stagate_label[i] = 'Layer_5'
    elif stagate_label[i] == '3':
        stagate_label[i] = 'Layer_4'
    elif stagate_label[i] == '1':
        stagate_label[i] = 'Layer_3'


specmix = pd.read_csv('SpecMix/151671/pred_label_151671.txt')
specmix_label = list(specmix['label SpiceMixPlus'])
specmix_label = [str(x) for x in specmix_label]
for i in range(len(specmix_label)):
    if specmix_label[i] == '3':
        specmix_label[i] = 'WM'
    elif specmix_label[i] == '0':
        specmix_label[i] = 'Layer_6'
    elif specmix_label[i] == '4':
        specmix_label[i] = 'Layer_5'
    elif specmix_label[i] == '2':
        specmix_label[i] = 'Layer_4'
    elif specmix_label[i] == '1':
        specmix_label[i] = 'Layer_3'

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


plot_color = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC']
# plot_color = ['#FABB6E', '#FC8002', '#ADDB88', '#369F2D', '#FAC7B3']
# plot_color = ['#6BB952', '#EC748B', '#C4A751', '#36ACA2', '#6EB1DE']
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


fig, axs = plt.subplots(3, 4, figsize=(13,10),constrained_layout=True)
sc.pl.spatial(adata, color=['Manual annotation'], ax=axs[0,0], show=False, spot_size=200, img_key=None, legend_fontsize=0, legend_loc='on data')
axs[0,0].spines['right'].set_visible(False) # 去掉边框
axs[0,0].spines['top'].set_visible(False)   # 去掉边框
axs[0,0].spines['left'].set_visible(False) # 去掉边框
axs[0,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('Manual annotation', fontsize=12)


sc.pl.spatial(adata, img_key=None, color=['sedr'], show=False, spot_size=200, ax=axs[0,1],
              legend_fontsize=0, legend_loc='on data')
axs[0,1].spines['right'].set_visible(False) # 去掉边框
axs[0,1].spines['top'].set_visible(False)   # 去掉边框
axs[0,1].spines['left'].set_visible(False) # 去掉边框
axs[0,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('SEDR  ARI:0.446', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['bayes'], show=False, spot_size=200, ax=axs[0,2],
              legend_fontsize=0, legend_loc='on data')
axs[0,2].spines['right'].set_visible(False) # 去掉边框
axs[0,2].spines['top'].set_visible(False)   # 去掉边框
axs[0,2].spines['left'].set_visible(False) # 去掉边框
axs[0,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,2].get_yaxis().set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].set_title('BayesSpace  ARI:0.733', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['spaceflow'], show=False, spot_size=200, ax=axs[0,3],
              legend_fontsize=0, legend_loc='on data')
axs[0,3].spines['right'].set_visible(False) # 去掉边框
axs[0,3].spines['top'].set_visible(False)   # 去掉边框
axs[0,3].spines['left'].set_visible(False) # 去掉边框
axs[0,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,3].get_yaxis().set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].set_title('SpaceFlow  ARI:0.327', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['deepst'], show=False, spot_size=200, ax=axs[1,0],
              legend_fontsize=0, legend_loc='on data')
axs[1,0].spines['right'].set_visible(False) # 去掉边框
axs[1,0].spines['top'].set_visible(False)   # 去掉边框
axs[1,0].spines['left'].set_visible(False) # 去掉边框
axs[1,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('DeepST  ARI:0.494', fontsize=12)


sc.pl.spatial(adata, img_key=None, color=['spagcn'], show=False, spot_size=200, ax=axs[1,1],
              legend_fontsize=0, legend_loc='on data')
axs[1,1].spines['right'].set_visible(False) # 去掉边框
axs[1,1].spines['top'].set_visible(False)   # 去掉边框
axs[1,1].spines['left'].set_visible(False) # 去掉边框
axs[1,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('SpaGCN  ARI:0.546', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['stlearn'], show=False, spot_size=200, ax=axs[1,2],
              legend_fontsize=0, legend_loc='on data')
axs[1,2].spines['right'].set_visible(False) # 去掉边框
axs[1,2].spines['top'].set_visible(False)   # 去掉边框
axs[1,2].spines['left'].set_visible(False) # 去掉边框
axs[1,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,2].get_yaxis().set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].set_title('stLearn  ARI:0.389', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['seurat'], show=False, spot_size=200, ax=axs[1,3],
              legend_fontsize=0, legend_loc='on data')
axs[1,3].spines['right'].set_visible(False) # 去掉边框
axs[1,3].spines['top'].set_visible(False)   # 去掉边框
axs[1,3].spines['left'].set_visible(False) # 去掉边框
axs[1,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,3].get_yaxis().set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].set_title('Seurat  ARI:0.273', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['stagate'], show=False, spot_size=200, ax=axs[2,0],
              legend_fontsize=0, legend_loc='on data')
axs[2,0].spines['right'].set_visible(False) # 去掉边框
axs[2,0].spines['top'].set_visible(False)   # 去掉边框
axs[2,0].spines['left'].set_visible(False) # 去掉边框
axs[2,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[2,0].get_yaxis().set_visible(False)
axs[2,0].get_xaxis().set_visible(False)
axs[2,0].set_title('STAGATE  ARI:0.359', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['scanpy'], show=False, spot_size=200, ax=axs[2,1],
              legend_fontsize=0, legend_loc='on data')
axs[2,1].spines['right'].set_visible(False) # 去掉边框
axs[2,1].spines['top'].set_visible(False)   # 去掉边框
axs[2,1].spines['left'].set_visible(False) # 去掉边框
axs[2,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[2,1].get_yaxis().set_visible(False)
axs[2,1].get_xaxis().set_visible(False)
axs[2,1].set_title('Scanpy  ARI:0.199', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['specmix'], show=False, spot_size=200, ax=axs[2,2],
              legend_fontsize=0, legend_loc='on data')
axs[2,2].spines['right'].set_visible(False) # 去掉边框
axs[2,2].spines['top'].set_visible(False)   # 去掉边框
axs[2,2].spines['left'].set_visible(False) # 去掉边框
axs[2,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2,2].get_yaxis().set_visible(False)
axs[2,2].get_xaxis().set_visible(False)
axs[2,2].set_title('SpecMix  ARI:0.458', fontsize=12)

sc.pl.spatial(adata, img_key=None, color=['PearlST'], show=False, spot_size=200, ax=axs[2,3]
                             )
axs[2,3].spines['right'].set_visible(False) # 去掉边框
axs[2,3].spines['top'].set_visible(False)   # 去掉边框
axs[2,3].spines['left'].set_visible(False) # 去掉边框
axs[2,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[2,3].get_yaxis().set_visible(False)
axs[2,3].get_xaxis().set_visible(False)
axs[2,3].set_title('PearlST  ARI:0.824', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/151671_domains.pdf', dpi=300)
# ari = adjusted_rand_score(true_label, stagate_label)



