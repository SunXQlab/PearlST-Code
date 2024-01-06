import scanpy as sc


def preprocessing_data(adata, n_top_genes):
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=5)
    sc.pp.normalize_total(adata, target_sum=1, exclude_highly_expressed=True, inplace=False)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    # adata = adata[:, adata.var.highly_variable]
    return adata