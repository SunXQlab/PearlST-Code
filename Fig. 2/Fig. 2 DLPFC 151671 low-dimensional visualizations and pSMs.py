import pandas as pd
import scanpy as sc
import numpy as np
from preprocessing import preprocessing_data
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
from utils import plot_pSM_multi, plot_pSM, prepare_figure_multi, pseudo_Spatiotemporal_Map


adata = sc.read_visium('151671')
adata = preprocessing_data(adata, n_top_genes=2000)

# load ground_truth
ground_truth = pd.read_csv('SpatialDE_clustering/cluster_labels_151671.csv', index_col=0)
ground_truth = list(ground_truth['ground_truth'])
adata.obs['true_label'] = pd.Categorical(ground_truth)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)
# sc.pl.umap(adata, color='true_label')
# sc.pl.spatial(adata, color='true_label')


adata_deepst = sc.read_h5ad('deepST/151671_deepst.h5ad')
deepst_label = list(adata_deepst.obs['DeepST_refine_domain'])
for i in range(len(deepst_label)):
    if deepst_label[i] == '4':
        deepst_label[i] = 'WM'
    elif deepst_label[i] == '3':
        deepst_label[i] = 'Layer_6'
    elif deepst_label[i] == '2':
        deepst_label[i] = 'Layer_5'
    elif deepst_label[i] == '1':
        deepst_label[i] = 'Layer_4'
    elif deepst_label[i] == '0':
        deepst_label[i] = 'Layer_3'
umap_deepst = anndata.AnnData(adata_deepst.obsm["DeepST_embed"])
umap_deepst.obs['deepst'] = deepst_label
sc.pp.neighbors(umap_deepst, n_neighbors=15)
sc.tl.umap(umap_deepst)
# sc.pl.umap(umap_deepst, color='DeepST_refine_domain')

sedr = pd.read_csv('SEDR/151671/metadata.tsv', index_col=0, sep='\t')
sedr_label = list(sedr['SEDR'])
sedr_label = [str(x) for x in sedr_label]
for i in range(len(sedr_label)):
    if sedr_label[i] == '4':
        sedr_label[i] = 'WM'
    elif sedr_label[i] == '3':
        sedr_label[i] = 'Layer_6'
    elif sedr_label[i] == '2':
        sedr_label[i] = 'Layer_5'
    elif sedr_label[i] == '1':
        sedr_label[i] = 'Layer_4'
    elif sedr_label[i] == '0':
        sedr_label[i] = 'Layer_3'
sedr_emb = np.load('SEDR/151671/SEDR_result.npz')
sedr_emb = sedr_emb['sedr_feat']
umap_sedr = anndata.AnnData(sedr_emb)
umap_sedr.obs['sedr'] = sedr_label
sc.pp.neighbors(umap_sedr, n_neighbors=15)
sc.tl.umap(umap_sedr)
# sc.pl.umap(umap_sedr, color='sedr')

spaceflow = sc.read_h5ad('spaceflow/151671_spaceflow.h5ad')
spaceflow_label = list(spaceflow.obs['leiden'])
spaceflow_label = [str(x) for x in spaceflow_label]
for i in range(len(spaceflow_label)):
    if spaceflow_label[i] == '4':
        spaceflow_label[i] = 'WM'
    elif spaceflow_label[i] == '3':
        spaceflow_label[i] = 'Layer_6'
    elif spaceflow_label[i] == '2':
        spaceflow_label[i] = 'Layer_5'
    elif spaceflow_label[i] == '1':
        spaceflow_label[i] = 'Layer_4'
    elif spaceflow_label[i] == '0':
        spaceflow_label[i] = 'Layer_3'
spaceflow_emb = spaceflow.obsm['embedding']
umap_spaceflow = anndata.AnnData(spaceflow_emb)
umap_spaceflow.obs['spaceflow'] = spaceflow_label
sc.pp.neighbors(umap_spaceflow, n_neighbors=15)
sc.tl.umap(umap_spaceflow)
# sc.pl.umap(umap_spaceflow, color='spaceflow')

stlearn = pd.read_csv('stlearn/151671/metadata.tsv', sep='\t')
stlearn_label = list(stlearn['X_pca_kmeans'])
stlearn_label = [str(x) for x in stlearn_label]
for i in range(len(stlearn_label)):
    if stlearn_label[i] == '4':
        stlearn_label[i] = 'WM'
    elif stlearn_label[i] == '3':
        stlearn_label[i] = 'Layer_6'
    elif stlearn_label[i] == '2':
        stlearn_label[i] = 'Layer_5'
    elif stlearn_label[i] == '1':
        stlearn_label[i] = 'Layer_4'
    elif stlearn_label[i] == '0':
        stlearn_label[i] = 'Layer_3'
stlearn_emb = pd.read_csv('stlearn/151671/PCs.tsv', sep='\t', index_col=0)
umap_stlearn = anndata.AnnData(stlearn_emb)
umap_stlearn.obs['stlearn'] = stlearn_label
sc.pp.neighbors(umap_stlearn, n_neighbors=15)
sc.tl.umap(umap_stlearn)
# sc.pl.umap(umap_stlearn, color='stlearn')

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
seurat_emb = pd.read_csv('Seurat/DLPFC/151671/Seurat/seurat.PCs.tsv', sep='\t')
umap_seurat = anndata.AnnData(seurat_emb)
umap_seurat.obs['seurat'] = seurat_label
sc.pp.neighbors(umap_seurat, n_neighbors=15)
sc.tl.umap(umap_seurat)


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
scanpy_emb = scanpy.obsm['X_pca']
umap_scanpy = anndata.AnnData(scanpy_emb)
umap_scanpy.obs['scanpy'] = scanpy_label
sc.pp.neighbors(umap_scanpy, n_neighbors=15)
sc.tl.umap(umap_scanpy)


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
stagate_emb = stagate.obsm['STAGATE']
umap_stagate = anndata.AnnData(stagate_emb)
umap_stagate.obs['stagate'] = stagate_label
sc.pp.neighbors(umap_stagate, n_neighbors=15)
sc.tl.umap(umap_stagate)


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
specmix_emb = np.load('SpecMix/151671/embedding151671.npy')
umap_specmix = anndata.AnnData(specmix_emb)
umap_specmix.obs['specmix'] = specmix_label
sc.pp.neighbors(umap_specmix, n_neighbors=15)
sc.tl.umap(umap_specmix)


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
warga_emb = warga.obsm['WARGA_embed']
umap_warga = anndata.AnnData(warga_emb)
umap_warga.obs['warga'] = warga_label
# umap_warga.obs['true_label'] = list(ground_truth['ground_truth'])
sc.pp.neighbors(umap_warga, n_neighbors=15)
sc.tl.umap(umap_warga)

ax = sc.pl.umap(umap_warga, color=['warga'], show=False,  size=30)
ax.spines['right'].set_visible(False) # 去掉边框
ax.spines['top'].set_visible(False)   # 去掉边框
ax.spines['left'].set_visible(False) # 去掉边框
ax.spines['bottom'].set_visible(False)   # 去掉边框
ax.get_yaxis().set_visible(False)
ax.get_xaxis().set_visible(False)
ax.set_title('PearlST', fontsize=0)
plt.savefig('embedding.pdf', dpi=300)


plot_color = ['#8ECFC9', '#FFBE7A', '#FA7F6F', '#82B0D2', '#BEB8DC']
umap_warga.uns['warga_colors'] = plot_color
umap_deepst.uns['deepst_colors'] = plot_color
umap_sedr.uns['sedr_colors'] = plot_color
umap_spaceflow.uns['spaceflow_colors'] = plot_color
umap_stlearn.uns['stlearn_colors'] = plot_color
umap_scanpy.uns['scanpy_colors'] = plot_color
umap_seurat.uns['seurat_colors'] = plot_color
umap_stagate.uns['stagate_colors'] = plot_color
umap_specmix.uns['specmix_colors'] = plot_color




fig, axs = plt.subplots(1, 4, figsize=(13,3),constrained_layout=True)
sc.pl.umap(umap_stagate, color=['stagate'],  ax=axs[0],
           show=False,  size=30, legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('STAGATE', fontsize=12)

sc.pl.umap(umap_spaceflow, color='spaceflow', ax=axs[1], show=False, size=30, legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('SpaceFlow', fontsize=12)

sc.pl.umap(umap_sedr, color='sedr', ax=axs[2], show=False, size=30, legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('SEDR', fontsize=12)

sc.pl.umap(umap_warga, color='warga', ax=axs[3], show=False, size=30, legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('PearlST', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/151671_embedding.pdf', dpi=300)



sedr_adata = adata.copy()
sedr_adata.obs['sedr'] = sedr_label
sedr_adata.obsm['SEDR_emb'] = sedr_emb
adata_seurat = adata.copy()
adata_seurat.obs['seurat'] = seurat_label
adata_seurat.obsm['seurat_emb'] = seurat_emb.values
adata_specmix = adata.copy()
adata_specmix.obs['specmix'] = specmix_label
adata_specmix.obsm['specmix_emb'] = specmix_emb

pseudo_Spatiotemporal_Map(sedr_adata, emb_name='SEDR_emb', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(warga, emb_name='WARGA_embed', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(stagate, emb_name='STAGATE', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(spaceflow, emb_name='embedding', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_deepst, emb_name='DeepST_embed', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(scanpy, emb_name='X_pca', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_seurat, emb_name='seurat_emb', n_neighbors=10, resolution=1.0)
pseudo_Spatiotemporal_Map(adata_specmix, emb_name='specmix_emb', n_neighbors=10, resolution=1.0)




plot_pSM_multi(adata=stagate, adata1=spaceflow, adata2=sedr_adata, adata3=warga, scatter_sz=30., rsz=3.,
             csz=13, wspace=.4, hspace=.5, left=0.008, right=0.98, bottom=0.1, top=0.9)


plt.savefig('C:/Users/haiyu/Desktop/151671_pSM_legend.pdf', dpi=300)
