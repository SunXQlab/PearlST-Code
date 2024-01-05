import anndata
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score
import matplotlib.pyplot as plt


adata = sc.read_h5ad('WARGA/stereoseq_chen_20_results.h5ad')
ari = adjusted_rand_score(adata.obs['annotation'], adata.obs['WARGA_refine_domain'])

labels = list(adata.obs['WARGA_refine_domain'])
for i in range(len(labels)):
    if labels[i] == '1':
        labels[i] = 'Connective tissue'
    elif labels[i] == '0':
        labels[i] = 'GI tract'
    elif labels[i] == '2':
        labels[i] = 'Epidermis'
    elif labels[i] == '3':
        labels[i] = 'Brain'
    elif labels[i] == '4':
        labels[i] = 'Heart'
    elif labels[i] == '5':
        labels[i] = 'Urogenital ridge'
    elif labels[i] == '6':
        labels[i] = 'Cartilage primordium'
    elif labels[i] == '7':
        labels[i] = 'Meninges'
    elif labels[i] == '8':
        labels[i] = 'Ovary'
    elif labels[i] == '9':
        labels[i] = 'Muscle'
    elif labels[i] == '10':
        labels[i] = 'Dorsal root ganglion'
    elif labels[i] == '11':
        labels[i] = 'Jaw and tooth'
    elif labels[i] == '12':
        labels[i] = 'Cavity'
    elif labels[i] == '13':
        labels[i] = 'Liver'
    elif labels[i] == '14':
        labels[i] = 'Meninges'
    elif labels[i] == '15':
        labels[i] = 'Muscle'
    elif labels[i] == '16':
        labels[i] = 'Cartilage primordium'
    elif labels[i] == '17':
        labels[i] = 'Kidney'
    elif labels[i] == '18':
        labels[i] = 'Brain'

adata.obs['WARGA_refine_domain'] = labels
# adata.uns['WARGA_refine_domain_colors'] = adata.uns['annotation_colors']
adata.uns['WARGA_refine_domain_colors'] = ['#ef833aff', '#3cb44bff', '#dfdce0ff', '#0bd3b1ff', '#b74c11ff',
                                           '#036df4ff', '#5c5ca6ff', '#d3245aff', '#f062f9ff', '#62cfe8ff',
                                           '#c923b1ff', '#dfca43ff', '#af1041ff', '#55afd9ff', '#9e7bffff']


fig, axs = plt.subplots(1, 2, figsize=(15,4),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial",  color='annotation', size=12, ax=axs[0])
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].invert_yaxis()
axs[0].set_title('Ground truth', fontsize=12)

sc.pl.embedding(adata, basis="spatial",  color='WARGA_refine_domain', size=12, ax=axs[1])
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].invert_yaxis()
axs[1].set_title('PearlST annotation', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/Stereoseq_chen_PearlST_annotation.pdf', dpi=300)


"""
sc.pl.embedding(adata, basis='spatial', groups='Muscle', color='annotation', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Liver', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Cartilage primordium', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Brain', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Cavity', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Dorsal root ganglion', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Epidermis', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Heart', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Meninges', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Muscle', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='Jaw and tooth', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis='spatial', groups='UNK', color='WARGA_refine_domain', size=13)
sc.pl.embedding(adata, basis="spatial", color='Map1b', size=15, legend_loc='on data', legend_fontsize=0)
sc.pl.embedding(adata, basis="spatial", color='Nnat', size=15, legend_loc='on data', legend_fontsize=0)
plt.savefig('stereo_70.jpg', dpi=300)
"""

fig, axs = plt.subplots(1, 8, figsize=(15,2.2),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Cartilage primordium',  color='annotation', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].invert_yaxis()
axs[0].set_title('Cartilage primordium', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Brain',  color='annotation', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].invert_yaxis()
axs[1].set_title('Brain', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Dorsal root ganglion',  color='annotation', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].invert_yaxis()
axs[2].set_title('Dorsal root ganglion', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Epidermis',  color='annotation', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].invert_yaxis()
axs[3].set_title('Epidermis', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Heart',  color='annotation', size=10, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].invert_yaxis()
axs[4].set_title('Heart', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Liver',  color='annotation', size=10, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].invert_yaxis()
axs[5].set_title('Liver', fontsize=12)


sc.pl.embedding(adata, basis="spatial", groups='Meninges',  color='annotation', size=10, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].invert_yaxis()
axs[6].set_title('Meninges', fontsize=12)

sc.pl.embedding(adata, basis="spatial", groups='Muscle',  color='annotation', size=10, ax=axs[7], legend_loc='on data', legend_fontsize=0)
axs[7].spines['right'].set_visible(False) # 去掉边框
axs[7].spines['top'].set_visible(False)   # 去掉边框
axs[7].spines['left'].set_visible(False) # 去掉边框
axs[7].spines['bottom'].set_visible(False)   # 去掉边框
axs[7].get_yaxis().set_visible(False)
axs[7].get_xaxis().set_visible(False)
axs[7].invert_yaxis()
axs[7].set_title('Muscle', fontsize=12)
plt.savefig('C:/Users/haiyu/Desktop/stereo_chen_groups_annotation.pdf', dpi=300)


fig, axs = plt.subplots(1, 8, figsize=(15,2.2), constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Cartilage primordium',  color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].invert_yaxis()
axs[0].set_title('Cartilage primordium', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Brain',  color='WARGA_refine_domain', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].invert_yaxis()
axs[1].set_title('Brain', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Dorsal root ganglion',  color='WARGA_refine_domain', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].invert_yaxis()
axs[2].set_title('Dorsal root ganglion', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Epidermis',  color='WARGA_refine_domain', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].invert_yaxis()
axs[3].set_title('Epidermis', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Heart',  color='WARGA_refine_domain', size=10, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].invert_yaxis()
axs[4].set_title('Heart', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Liver',  color='WARGA_refine_domain', size=10, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].invert_yaxis()
axs[5].set_title('Liver', fontsize=0)


sc.pl.embedding(adata, basis="spatial", groups='Meninges',  color='WARGA_refine_domain', size=10, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].invert_yaxis()
axs[6].set_title('Meninges', fontsize=0)

sc.pl.embedding(adata, basis="spatial", groups='Muscle',  color='WARGA_refine_domain', size=10, ax=axs[7], legend_loc='on data', legend_fontsize=0)
axs[7].spines['right'].set_visible(False) # 去掉边框
axs[7].spines['top'].set_visible(False)   # 去掉边框
axs[7].spines['left'].set_visible(False) # 去掉边框
axs[7].spines['bottom'].set_visible(False)   # 去掉边框
axs[7].get_yaxis().set_visible(False)
axs[7].get_xaxis().set_visible(False)
axs[7].invert_yaxis()
axs[7].set_title('Muscle', fontsize=0)

plt.savefig('C:/Users/haiyu/Desktop/stereo_chen_groups_PearlST.pdf', dpi=300)



fig, axs = plt.subplots(1, 8, figsize=(15,2.0),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", color='Meox1', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].invert_yaxis()
axs[0].set_title('Meox1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Tubb2b', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].invert_yaxis()
axs[1].set_title('Tubb2b', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Ntrk1', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].invert_yaxis()
axs[2].set_title('Ntrk1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Krt5', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].invert_yaxis()
axs[3].set_title('Krt5', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Myh7', size=10, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].invert_yaxis()
axs[4].set_title('Myh7', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Alb', size=10, ax=axs[5], legend_loc='on data', legend_fontsize=0)
axs[5].spines['right'].set_visible(False) # 去掉边框
axs[5].spines['top'].set_visible(False)   # 去掉边框
axs[5].spines['left'].set_visible(False) # 去掉边框
axs[5].spines['bottom'].set_visible(False)   # 去掉边框
axs[5].get_yaxis().set_visible(False)
axs[5].get_xaxis().set_visible(False)
axs[5].invert_yaxis()
axs[5].set_title('Alb', fontsize=12)


sc.pl.embedding(adata, basis="spatial", color='Map1b', size=10, ax=axs[6], legend_loc='on data', legend_fontsize=0)
axs[6].spines['right'].set_visible(False) # 去掉边框
axs[6].spines['top'].set_visible(False)   # 去掉边框
axs[6].spines['left'].set_visible(False) # 去掉边框
axs[6].spines['bottom'].set_visible(False)   # 去掉边框
axs[6].get_yaxis().set_visible(False)
axs[6].get_xaxis().set_visible(False)
axs[6].invert_yaxis()
axs[6].set_title('Map1b', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Myog', size=10, ax=axs[7], legend_loc='on data', legend_fontsize=0)
axs[7].spines['right'].set_visible(False) # 去掉边框
axs[7].spines['top'].set_visible(False)   # 去掉边框
axs[7].spines['left'].set_visible(False) # 去掉边框
axs[7].spines['bottom'].set_visible(False)   # 去掉边框
axs[7].get_yaxis().set_visible(False)
axs[7].get_xaxis().set_visible(False)
axs[7].invert_yaxis()
axs[7].set_title('Myog', fontsize=12)

plt.savefig('C:/Users/haiyu/Desktop/stereo_chen_groups_markers.pdf', dpi=300)



fig, axs = plt.subplots(1, 4, figsize=(13,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Cartilage primordium', color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Cartilage primordium', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Meox1', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Meox1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Pax1', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Pax1', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Pax9', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Pax9', fontsize=12)
plt.savefig('Stereoseq_chen_Cartilage primordium_markers.pdf', dpi=300)


fig, axs = plt.subplots(1, 3, figsize=(9.5,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Dorsal root ganglion', color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Dorsal root ganglion', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Prdm12', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Prdm12', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Ntrk1', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Ntrk1', fontsize=12)
plt.savefig('Stereo_chen_Dorsal root ganglion_markers.pdf', dpi=300)



fig, axs = plt.subplots(1, 4, figsize=(13,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Epidermis', color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Epidermis', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Krt5', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Krt5', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Krt14', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Krt14', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Krt15', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Krt15', fontsize=12)
plt.savefig('Stereo_chen_Epidermis primordium_markers.pdf', dpi=300)



fig, axs = plt.subplots(1, 5, figsize=(16,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Heart', color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Heart', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Nppa', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Nppa', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Myl7', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Myl7', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Myl2', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Myl2', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Myh7', size=10, ax=axs[4], legend_loc='on data', legend_fontsize=0)
axs[4].spines['right'].set_visible(False) # 去掉边框
axs[4].spines['top'].set_visible(False)   # 去掉边框
axs[4].spines['left'].set_visible(False) # 去掉边框
axs[4].spines['bottom'].set_visible(False)   # 去掉边框
axs[4].get_yaxis().set_visible(False)
axs[4].get_xaxis().set_visible(False)
axs[4].set_title('Myh7', fontsize=12)
plt.savefig('Stereo_chen_Heart_markers.pdf', dpi=300)



fig, axs = plt.subplots(1, 4, figsize=(13,3),constrained_layout=True)
sc.pl.embedding(adata, basis="spatial", groups='Liver', color='WARGA_refine_domain', size=10, ax=axs[0], legend_loc='on data', legend_fontsize=0)
axs[0].spines['right'].set_visible(False) # 去掉边框
axs[0].spines['top'].set_visible(False)   # 去掉边框
axs[0].spines['left'].set_visible(False) # 去掉边框
axs[0].spines['bottom'].set_visible(False)   # 去掉边框
axs[0].get_yaxis().set_visible(False)
axs[0].get_xaxis().set_visible(False)
axs[0].set_title('Liver', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Afp', size=10, ax=axs[1], legend_loc='on data', legend_fontsize=0)
axs[1].spines['right'].set_visible(False) # 去掉边框
axs[1].spines['top'].set_visible(False)   # 去掉边框
axs[1].spines['left'].set_visible(False) # 去掉边框
axs[1].spines['bottom'].set_visible(False)   # 去掉边框
axs[1].get_yaxis().set_visible(False)
axs[1].get_xaxis().set_visible(False)
axs[1].set_title('Afp', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Alb', size=10, ax=axs[2], legend_loc='on data', legend_fontsize=0)
axs[2].spines['right'].set_visible(False) # 去掉边框
axs[2].spines['top'].set_visible(False)   # 去掉边框
axs[2].spines['left'].set_visible(False) # 去掉边框
axs[2].spines['bottom'].set_visible(False)   # 去掉边框
axs[2].get_yaxis().set_visible(False)
axs[2].get_xaxis().set_visible(False)
axs[2].set_title('Alb', fontsize=12)

sc.pl.embedding(adata, basis="spatial", color='Fgb', size=10, ax=axs[3], legend_loc='on data', legend_fontsize=0)
axs[3].spines['right'].set_visible(False) # 去掉边框
axs[3].spines['top'].set_visible(False)   # 去掉边框
axs[3].spines['left'].set_visible(False) # 去掉边框
axs[3].spines['bottom'].set_visible(False)   # 去掉边框
axs[3].get_yaxis().set_visible(False)
axs[3].get_xaxis().set_visible(False)
axs[3].set_title('Fgb', fontsize=12)
plt.savefig('Stereo_chen_Liver_markers.pdf', dpi=300)



