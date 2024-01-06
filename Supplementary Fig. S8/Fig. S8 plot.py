import scanpy as sc
import anndata
import pandas as pd
from utils import pseudo_Spatiotemporal_Map, plot_pSM, plot_pSM_multi
import matplotlib.pyplot as plt


adata1 = sc.read_h5ad('./WARGA/DLPFC dataset/augment test/151507.h5ad')     # no augment
adata1_aug = sc.read_h5ad('./WARGA/DLPFC dataset/151507/PearlST_results.h5ad')   # augment

adata2 = sc.read_h5ad('./WARGA/DLPFC dataset/augment test/151509.h5ad')     # no augment
adata2_aug = sc.read_h5ad('./WARGA/DLPFC dataset/151509/PearlST_results.h5ad')   # augment

adata3 = sc.read_h5ad('./WARGA/DLPFC dataset/augment test/151670.h5ad')     # no augment
adata3_aug = sc.read_h5ad('./WARGA/DLPFC dataset/151670/PearlST_results.h5ad')   # augment

adata4 = sc.read_h5ad('./WARGA/DLPFC dataset/augment test/151672.h5ad')     # no augment
adata4_aug = sc.read_h5ad('./WARGA/DLPFC dataset/151672/PearlST_results.h5ad')   # augment

adata5 = sc.read_h5ad('./WARGA/DLPFC dataset/augment test/151675.h5ad')     # no augment
adata5_aug = sc.read_h5ad('./WARGA/DLPFC dataset/151675/PearlST_results.h5ad')   # augment


# low-dimensional visualization
umap_adata1 = anndata.AnnData(adata1.obsm['WARGA_embed'])
umap_adata1.obs['pred_label'] = list(adata1.obs['pred_label'])
sc.pp.neighbors(umap_adata1, n_neighbors=10)
sc.tl.umap(umap_adata1)
sc.tl.paga(umap_adata1, groups='pred_label')
# sc.pl.umap(umap_adata1, color='pred_label')

umap_adata1_aug = anndata.AnnData(adata1_aug.obsm['WARGA_embed'])
umap_adata1_aug.obs['pred_label'] = list(adata1_aug.obs['pred_label'])
sc.pp.neighbors(umap_adata1_aug, n_neighbors=10)
sc.tl.umap(umap_adata1_aug)
# sc.pl.umap(umap_adata1_aug, color='pred_label')
sc.tl.paga(umap_adata1_aug, groups='pred_label')


umap_adata2 = anndata.AnnData(adata2.obsm['WARGA_embed'])
umap_adata2.obs['pred_label'] = list(adata2.obs['pred_label'])
sc.pp.neighbors(umap_adata2, n_neighbors=10)
sc.tl.umap(umap_adata2)
# sc.pl.umap(umap_adata2, color='pred_label')
sc.tl.paga(umap_adata2, groups='pred_label')

umap_adata2_aug = anndata.AnnData(adata2_aug.obsm['WARGA_embed'])
umap_adata2_aug.obs['pred_label'] = list(adata2_aug.obs['pred_label'])
sc.pp.neighbors(umap_adata2_aug, n_neighbors=10)
sc.tl.umap(umap_adata2_aug)
# sc.pl.umap(umap_adata2_aug, color='pred_label')
sc.tl.paga(umap_adata2_aug, groups='pred_label')


umap_adata3 = anndata.AnnData(adata3.obsm['WARGA_embed'])
umap_adata3.obs['pred_label'] = list(adata3.obs['pred_label'])
sc.pp.neighbors(umap_adata3, n_neighbors=10)
sc.tl.umap(umap_adata3)
# sc.pl.umap(umap_adata3, color='pred_label')
sc.tl.paga(umap_adata3, groups='pred_label')

umap_adata3_aug = anndata.AnnData(adata3_aug.obsm['WARGA_embed'])
umap_adata3_aug.obs['pred_label'] = list(adata3_aug.obs['pred_label'])
sc.pp.neighbors(umap_adata3_aug, n_neighbors=10)
sc.tl.umap(umap_adata3_aug)
# sc.pl.umap(umap_adata3_aug, color='pred_label')
sc.tl.paga(umap_adata3_aug, groups='pred_label')


umap_adata4 = anndata.AnnData(adata4.obsm['WARGA_embed'])
umap_adata4.obs['pred_label'] = list(adata4.obs['pred_label'])
sc.pp.neighbors(umap_adata4, n_neighbors=10)
sc.tl.umap(umap_adata4)
# sc.pl.umap(umap_adata4, color='pred_label')
sc.tl.paga(umap_adata4, groups='pred_label')

umap_adata4_aug = anndata.AnnData(adata4_aug.obsm['WARGA_embed'])
umap_adata4_aug.obs['pred_label'] = list(adata4_aug.obs['pred_label'])
sc.pp.neighbors(umap_adata4_aug, n_neighbors=10)
sc.tl.umap(umap_adata4_aug)
# sc.pl.umap(umap_adata4_aug, color='pred_label')
sc.tl.paga(umap_adata4_aug, groups='pred_label')


umap_adata5 = anndata.AnnData(adata5.obsm['WARGA_embed'])
umap_adata5.obs['pred_label'] = list(adata5.obs['pred_label'])
sc.pp.neighbors(umap_adata5, n_neighbors=10)
sc.tl.umap(umap_adata5)
# sc.pl.umap(umap_adata5, color='pred_label')
sc.tl.paga(umap_adata5, groups='pred_label')

umap_adata5_aug = anndata.AnnData(adata5_aug.obsm['WARGA_embed'])
umap_adata5_aug.obs['pred_label'] = list(adata5_aug.obs['pred_label'])
sc.pp.neighbors(umap_adata5_aug, n_neighbors=10)
sc.tl.umap(umap_adata5_aug)
# sc.pl.umap(umap_adata5_aug, color='pred_label')
sc.tl.paga(umap_adata5_aug, groups='pred_label')

# visualizations
fig, axs = plt.subplots(2, 5, figsize=(13,5),constrained_layout=True)
sc.pl.umap(umap_adata1, color='pred_label', size=12, ax=axs[0,0], legend_loc='on data', legend_fontsize=0)
axs[0,0].spines['right'].set_visible(False) # 去掉边框
axs[0,0].spines['top'].set_visible(False)   # 去掉边框
axs[0,0].spines['left'].set_visible(False) # 去掉边框
axs[0,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata2, color='pred_label', size=12, ax=axs[0,1], legend_loc='on data', legend_fontsize=0)
axs[0,1].spines['right'].set_visible(False) # 去掉边框
axs[0,1].spines['top'].set_visible(False)   # 去掉边框
axs[0,1].spines['left'].set_visible(False) # 去掉边框
axs[0,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata3, color='pred_label', size=12, ax=axs[0,2], legend_loc='on data', legend_fontsize=0)
axs[0,2].spines['right'].set_visible(False) # 去掉边框
axs[0,2].spines['top'].set_visible(False)   # 去掉边框
axs[0,2].spines['left'].set_visible(False) # 去掉边框
axs[0,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,2].get_yaxis().set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata4, color='pred_label', size=12, ax=axs[0,3], legend_loc='on data', legend_fontsize=0)
axs[0,3].spines['right'].set_visible(False) # 去掉边框
axs[0,3].spines['top'].set_visible(False)   # 去掉边框
axs[0,3].spines['left'].set_visible(False) # 去掉边框
axs[0,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,3].get_yaxis().set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata5, color='pred_label', size=12, ax=axs[0,4], legend_loc='on data', legend_fontsize=0)
axs[0,4].spines['right'].set_visible(False) # 去掉边框
axs[0,4].spines['top'].set_visible(False)   # 去掉边框
axs[0,4].spines['left'].set_visible(False) # 去掉边框
axs[0,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,4].get_yaxis().set_visible(False)
axs[0,4].get_xaxis().set_visible(False)
axs[0,4].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata1_aug, color='pred_label', size=12, ax=axs[1,0], legend_loc='on data', legend_fontsize=0)
axs[1,0].spines['right'].set_visible(False) # 去掉边框
axs[1,0].spines['top'].set_visible(False)   # 去掉边框
axs[1,0].spines['left'].set_visible(False) # 去掉边框
axs[1,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata2_aug, color='pred_label', size=12, ax=axs[1,1], legend_loc='on data', legend_fontsize=0)
axs[1,1].spines['right'].set_visible(False) # 去掉边框
axs[1,1].spines['top'].set_visible(False)   # 去掉边框
axs[1,1].spines['left'].set_visible(False) # 去掉边框
axs[1,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata3_aug, color='pred_label', size=12, ax=axs[1,2], legend_loc='on data', legend_fontsize=0)
axs[1,2].spines['right'].set_visible(False) # 去掉边框
axs[1,2].spines['top'].set_visible(False)   # 去掉边框
axs[1,2].spines['left'].set_visible(False) # 去掉边框
axs[1,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,2].get_yaxis().set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata4_aug, color='pred_label', size=12, ax=axs[1,3], legend_loc='on data', legend_fontsize=0)
axs[1,3].spines['right'].set_visible(False) # 去掉边框
axs[1,3].spines['top'].set_visible(False)   # 去掉边框
axs[1,3].spines['left'].set_visible(False) # 去掉边框
axs[1,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,3].get_yaxis().set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].set_title('pred_label', fontsize=0)

sc.pl.umap(umap_adata5_aug, color='pred_label', size=12, ax=axs[1,4], legend_loc='on data', legend_fontsize=0)
axs[1,4].spines['right'].set_visible(False) # 去掉边框
axs[1,4].spines['top'].set_visible(False)   # 去掉边框
axs[1,4].spines['left'].set_visible(False) # 去掉边框
axs[1,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,4].get_yaxis().set_visible(False)
axs[1,4].get_xaxis().set_visible(False)
axs[1,4].set_title('pred_label', fontsize=0)
plt.savefig('./augment test/low-dimensional viaualizations.pdf', dpi=300)


# paga visualizations
fig, axs = plt.subplots(2, 5, figsize=(13,5),constrained_layout=True)
sc.pl.paga(umap_adata1, color='pred_label', ax=axs[0,0], show=False, fontsize=10)
axs[0,0].spines['right'].set_visible(False) # 去掉边框
axs[0,0].spines['top'].set_visible(False)   # 去掉边框
axs[0,0].spines['left'].set_visible(False) # 去掉边框
axs[0,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,0].get_yaxis().set_visible(False)
axs[0,0].get_xaxis().set_visible(False)
axs[0,0].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata2, color='pred_label', ax=axs[0,1], show=False, fontsize=10)
axs[0,1].spines['right'].set_visible(False) # 去掉边框
axs[0,1].spines['top'].set_visible(False)   # 去掉边框
axs[0,1].spines['left'].set_visible(False) # 去掉边框
axs[0,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,1].get_yaxis().set_visible(False)
axs[0,1].get_xaxis().set_visible(False)
axs[0,1].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata3, color='pred_label', ax=axs[0,2], show=False, fontsize=10)
axs[0,2].spines['right'].set_visible(False) # 去掉边框
axs[0,2].spines['top'].set_visible(False)   # 去掉边框
axs[0,2].spines['left'].set_visible(False) # 去掉边框
axs[0,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,2].get_yaxis().set_visible(False)
axs[0,2].get_xaxis().set_visible(False)
axs[0,2].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata4, color='pred_label', ax=axs[0,3], show=False, fontsize=10)
axs[0,3].spines['right'].set_visible(False) # 去掉边框
axs[0,3].spines['top'].set_visible(False)   # 去掉边框
axs[0,3].spines['left'].set_visible(False) # 去掉边框
axs[0,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,3].get_yaxis().set_visible(False)
axs[0,3].get_xaxis().set_visible(False)
axs[0,3].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata5, color='pred_label', ax=axs[0,4], show=False, fontsize=10)
axs[0,4].spines['right'].set_visible(False) # 去掉边框
axs[0,4].spines['top'].set_visible(False)   # 去掉边框
axs[0,4].spines['left'].set_visible(False) # 去掉边框
axs[0,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[0,4].get_yaxis().set_visible(False)
axs[0,4].get_xaxis().set_visible(False)
axs[0,4].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata1_aug, color='pred_label', ax=axs[1,0], show=False, fontsize=10)
axs[1,0].spines['right'].set_visible(False) # 去掉边框
axs[1,0].spines['top'].set_visible(False)   # 去掉边框
axs[1,0].spines['left'].set_visible(False) # 去掉边框
axs[1,0].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,0].get_yaxis().set_visible(False)
axs[1,0].get_xaxis().set_visible(False)
axs[1,0].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata2_aug, color='pred_label', ax=axs[1,1], show=False, fontsize=10)
axs[1,1].spines['right'].set_visible(False) # 去掉边框
axs[1,1].spines['top'].set_visible(False)   # 去掉边框
axs[1,1].spines['left'].set_visible(False) # 去掉边框
axs[1,1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,1].get_yaxis().set_visible(False)
axs[1,1].get_xaxis().set_visible(False)
axs[1,1].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata3_aug, color='pred_label', ax=axs[1,2], show=False, fontsize=10)
axs[1,2].spines['right'].set_visible(False) # 去掉边框
axs[1,2].spines['top'].set_visible(False)   # 去掉边框
axs[1,2].spines['left'].set_visible(False) # 去掉边框
axs[1,2].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,2].get_yaxis().set_visible(False)
axs[1,2].get_xaxis().set_visible(False)
axs[1,2].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata4_aug, color='pred_label', ax=axs[1,3], show=False, fontsize=10)
axs[1,3].spines['right'].set_visible(False) # 去掉边框
axs[1,3].spines['top'].set_visible(False)   # 去掉边框
axs[1,3].spines['left'].set_visible(False) # 去掉边框
axs[1,3].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,3].get_yaxis().set_visible(False)
axs[1,3].get_xaxis().set_visible(False)
axs[1,3].set_title('pred_label', fontsize=0)

sc.pl.paga(umap_adata5_aug, color='pred_label', ax=axs[1,4], show=False, fontsize=10)
axs[1,4].spines['right'].set_visible(False) # 去掉边框
axs[1,4].spines['top'].set_visible(False)   # 去掉边框
axs[1,4].spines['left'].set_visible(False) # 去掉边框
axs[1,4].spines['bottom'].set_visible(False)   # 去掉边框
axs[1,4].get_yaxis().set_visible(False)
axs[1,4].get_xaxis().set_visible(False)
axs[1,4].set_title('pred_label', fontsize=0)
plt.savefig('./augment test/paga visualizations.pdf', dpi=300)



# pSM visualization
pseudo_Spatiotemporal_Map(adata1, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata1_aug, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata2, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata2_aug, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata3, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata3_aug, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata4, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata4_aug, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata5, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)
pseudo_Spatiotemporal_Map(adata5_aug, emb_name='WARGA_embed', n_neighbors=20, resolution=1.0)

plot_pSM_multi(adata=adata1, adata1=adata2, adata2=adata3, adata3=adata4, adata4=adata5, scatter_sz=30., rsz=3.,
             csz=13, wspace=.4, hspace=.5, left=0.008, right=0.98, bottom=0.1, top=0.9)
plt.savefig('./augment test/pSM no-augment.pdf', dpi=300)

plot_pSM_multi(adata=adata1_aug, adata1=adata2_aug, adata2=adata3_aug, adata3=adata4_aug, adata4=adata5_aug, scatter_sz=30., rsz=3.,
             csz=13, wspace=.4, hspace=.5, left=0.008, right=0.98, bottom=0.1, top=0.9)
plt.savefig('./augment test/pSM augment.pdf', dpi=300)

plot_pSM(adata1_aug, scatter_sz=20., rsz=4.,
             csz=4.5, wspace=.4, hspace=.5, left=0.125, right=0.95, bottom=0.1, top=0.9)
plt.savefig('./augment test/pSM legend.pdf', dpi=300)


