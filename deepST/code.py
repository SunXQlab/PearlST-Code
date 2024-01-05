import scanpy as sc


adata = sc.read_h5ad('./deepST/Human breast cancer_deepst_results.h5ad')

del adata.obsm['adjacent_data']
del adata.obsm['augment_gene_data']
del adata.obsm['image_feat']
del adata.obsm['image_feat_pca']
del adata.obsm['weights_matrix_all']

adata.write_h5ad('./deepST/Human breast cancer_deepst_results.h5ad', compression='gzip')