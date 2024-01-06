import scanpy as sc
import numpy as np
import h5py
import pandas as pd
import matplotlib.pyplot as plt
import plotting
import seaborn as sns
from WARGA.model import GCNModelAE
from WARGA.preprocessing import preprocessing_data, graph_alpha
import torch
import scipy.sparse as sp
from WARGA.utils import preprocess_graph
import scipy.stats as stats


# load clustering results obtained by PearlST
adata_pearlST = sc.read_h5ad('./WARGA/DLPFC dataset/Human breast cancer/PearlST_annotation.h5ad')

# load Ligand-Receptor Database
LR_DB = h5py.File('./WARGA/DLPFC dataset/human breast cancer/LigRec_DB.h5')
Ligand = list(LR_DB['Lig'])
Ligand = [x.astype(np.str_) for x in Ligand]
Receptor = list(LR_DB['Rec'])
Receptor = [x.astype(np.str_) for x in Receptor]

# read raw data
# adata = sc.read_visium('./WARGA/DLPFC dataset/human breast cancer')
adata = adata_pearlST
adata.X = adata.obsm['augment_data']


# adata.obs['annotation'] = adata_pearlST.obs['annotation']

"""
Finding all L-R pairs expressed in gene expression matrix
"""
def find_all_LR_pair_expressed(adata, Ligand, Receptor):
    LR_pairs = []
    for i in range(len(Ligand)):
        if Ligand[i] in adata.var_names and Receptor[i] in adata.var_names:
            LR_pairs.append([Ligand[i], Receptor[i]])
    return LR_pairs


def LR_pair_activated(adata_1, adata_2, LR_pair):
    """
    param adata_1: ligand cells expression
    param adata_2: receptor cells expression
    LR_pair: selected ligand-receptor pair
    return: the signal strength received by receptor cells
    """

    LR_score_pair = np.zeros(shape=adata_2.shape[0])
    for j in range(adata_2.shape[0]):
        # compute distances between each ligand cell and jth receptor cell
        dis = np.linalg.norm(adata_1.obsm['spatial'] - adata_2.obsm['spatial'][j], axis=1)
        dis = dis[:, np.newaxis]
        dis_rec = np.reciprocal(dis)
        LR_j = np.sum(adata_1[:, LR_pair[0]].X.toarray() * dis_rec * adata_2[:, LR_pair[1]].X[j].toarray())
        LR_score_pair[j] = LR_j
    return LR_score_pair



def compute_LR_score(adata, LR_pairs, lig_cell, rec_cell):
    """
    param adata: whole gene expression and spatial locations
    param LR_pairs: selected all LR-pairs from prior L-R database
    param lig_cell: the cell type or region of selected ligand cells
    param rec_cell: the cell type or region of selected receptor cells
    return: LR-score matrix, which represents the signal strength received by receptor cells on all LR-pairs
    """

    adata_lig = adata[adata.obs['annotation'] == lig_cell]
    adata_rec = adata[adata.obs['annotation'] == rec_cell]

    LR_score = np.zeros(shape=(adata_rec.shape[0], len(LR_pairs)))
    for i in range(len(LR_pairs)):
        LR_score_pair_i = LR_pair_activated(adata_lig, adata_rec, LR_pairs[i])
        LR_score[:, i] = LR_score_pair_i
        print("Ligand-Receptor pair " + str(i) + ":" + LR_pairs[i][0] + "-" + LR_pairs[i][1])
    return LR_score


LR_pairs = find_all_LR_pair_expressed(adata, Ligand, Receptor)
LR_score = compute_LR_score(adata, LR_pairs, lig_cell='IDC_4', rec_cell='Tumor_edge_1')
lr_pair_name = [x[0] + '-' + x[1] for x in LR_pairs]
LR_score = pd.DataFrame(LR_score, columns=lr_pair_name)
LR_score.to_csv('./CCI/IDC_4-Tumor_edge_1 LR_score.csv')


LR_score = pd.read_csv('./CCI/IDC_4-Tumor_edge_1 LR_score.csv', index_col=0)
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

# preparing training data
X = LR_score
y = adata[adata.obs['annotation'] == 'Tumor_edge_1'].obsm['WARGA_embed']


# construct random forest regression model and set multi-cores to perform parallel computation
rf = RandomForestRegressor(n_estimators=400, n_jobs=6, random_state=0)
# fitting model using training data
rf.fit(X, y)
# test on testing data
y_pred = rf.predict(X)
# evaluate model performance
mse = mean_squared_error(y, y_pred)
print("均方误差(MSE):", mse)

plt.figure()
plt.scatter(y_pred, y, alpha=0.2)
plotting.lineOfIdentity()
plotting.addCorrelation(y_pred, y)
plt.xlabel('Fit')
plt.ylabel('Data')
plt.gca().axis('equal')
plt.gca().set_xticks([0, 0.5, 1])
plt.gca().set_yticks([0, 0.5, 1])
plt.savefig('CCI/IDC_4_Tumor_edge_1_RF_fitting.pdf', dpi=300)

outNameGene = ['embedding_' + str(i) for i in range(32)]
plt.rcParams["figure.figsize"] = (10,9)
plt.figure()
rank = plotting.compareAllTFs(y_pred, y, outNameGene)
plt.savefig('CCI/IDC_4_Tumor_edge_1_RF_embeddings.pdf', dpi=300)

"""
Find low-dimensional representations that have significant differences
"""
def find_variable_embeddings(adata, top_k):
    data = adata.obsm['WARGA_embed']
    new_dict = {}
    for i in range(data.shape[1]):
        key_i = "embedding_{}".format(i)
        feature_i = data[:, i]
        new_dict[key_i] = {"features": feature_i, "p-times": []}

    # wilcoxon test
    keys = list(new_dict.keys())
    for i in range(len(new_dict)):
        for j in range(i+1, len(new_dict)):
            statistic, p_value = stats.wilcoxon(new_dict[keys[i]]['features'], new_dict[keys[j]]['features'])
            if p_value < 0.05:
                new_dict[keys[i]]['p-times'].append(str(i) + "-" + str(j))
                new_dict[keys[j]]['p-times'].append(str(i) + "-" + str(j))

    # for i in range(len(new_dict)):
    #     print(len(new_dict[keys[i]]['p-times']))

    sorted_dict = sorted(new_dict.items(), key=lambda x: len(x[1]['p-times']), reverse=True)

    # for i in range(len(sorted_dict)):
    #     print(len(sorted_dict[i][1]['p-times']))

    top_k_keys = [key for key, value in sorted_dict[:top_k]]
    return top_k_keys


# get variable importance using "feature importance" from random forest regression
def compute_feature_importance(emb_target, outNameGene):
    target_index = outNameGene.index(emb_target)  # target index
    target_estimator = rf.estimators_[target_index]
    # get feature importance of specific target
    target_importances = target_estimator.feature_importances_

    feature_names = X.columns     # LR-score
    importance_df = pd.DataFrame({'features': feature_names, 'importance': target_importances})
    importance_df = importance_df.sort_values(by='importance', ascending=False)
    return importance_df, target_index


"""
Find genes that have major contributions for top-k embeddings having significant differences
"""
def get_top_important_genes_for_target_emb(adata, file_path, target_index, top_k):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = GCNModelAE(2000, 256, 32, 0.1).to(device)
    model.load_state_dict(torch.load(file_path))
    model.eval()

    features = adata.obsm['augment_data']
    n_nodes, feat_dim = features.shape
    features = torch.tensor(features).to(device)
    adj = graph_alpha(adata_pearlST.obsm['spatial'], n_neighbors=10)
    adj = adj.astype(np.float32)
    # Store original adjacency matrix (without diagonal entries) for later
    adj_orig = adj
    adj_orig = adj_orig - sp.dia_matrix((adj_orig.diagonal()[np.newaxis, :], [0]), shape=adj_orig.shape)
    adj_orig.eliminate_zeros()
    # Some preprocessing
    adj_norm = preprocess_graph(adj_orig)
    adj_norm = adj_norm.to(device)
    _, baseline_predictions = model.forward(x=features, adj=adj_norm)
    baseline_predictions = baseline_predictions.cpu().detach().numpy()

    target_emb_importance = []
    for i in range(features.shape[1]):
        X_perturbed = features.clone()
        X_perturbed[:, i] = 0
        _, perturbed_predictions = model.forward(x=X_perturbed, adj=adj_norm)
        perturbed_predictions = perturbed_predictions.cpu().detach().numpy()
        feature_importance = np.mean(np.abs(perturbed_predictions[:, target_index] - baseline_predictions[:, target_index]))
        target_emb_importance.append(feature_importance)
    target_emb_importance_sort = sorted(target_emb_importance, reverse=True)
    top_important_index = [target_emb_importance.index(x) for x in target_emb_importance_sort[:top_k]]
    top_important_genes = adata.var_names[top_important_index]
    return top_important_genes, target_emb_importance_sort[:top_k]

top_5_emb = find_variable_embeddings(adata, top_k=2)


def construct_multi_network(top_emb):
    net = {}
    for emb in top_emb:
        importance_df, target_index = compute_feature_importance(emb, outNameGene)
        plt.figure(figsize=(10, 8))
        sns.barplot(x='features', y='importance', data=importance_df[:2])
        plt.xticks(rotation=30)
        plt.xlabel("LR-scores")
        plt.ylabel("Importances")
        plt.title("feature importance about embedding_{}".format(target_index))
        plt.show()
        # plt.savefig('feature importance about embedding_{}.jpg'.format(target_index), dpi=300)
        top_important_genes, top_importances = get_top_important_genes_for_target_emb(adata,
                                                            file_path='DLPFC dataset/Human breast cancer/model.pt',
                                                             target_index=target_index, top_k=3)
        net[emb] = {"LR_pairs": importance_df[:2], "target_genes": top_important_genes, "importance": top_importances}
    return net


net = construct_multi_network(top_emb=top_5_emb)


import pickle
file_path = './CCI/IDC_4-Tumor_edge_1 net_workflow.pkl'
# 使用pickle保存字典到文件
with open(file_path, 'wb') as file:
    pickle.dump(net, file)

plt.savefig('CCI/IDC_4_Tumor_edge_1_LR_importance_embedding_7.pdf', dpi=300)
plt.savefig('CCI/IDC_4_Tumor_edge_1_LR_importance_embedding_12.pdf', dpi=300)
plt.savefig('CCI/IDC_4_Tumor_edge_1_LR_importance_embedding_9.pdf', dpi=300)
plt.savefig('CCI/IDC_4_Tumor_edge_1_LR_importance_embedding_5.pdf', dpi=300)
plt.savefig('CCI/IDC_4_Tumor_edge_1_LR_importance_embedding_6.pdf', dpi=300)


"""
# Gene enrichment analysis
from gprofiler import GProfiler
gp=GProfiler(user_agent='ExampleTool', return_dataframe=True)
gene_list = ["ENSG00000157764", "ENSG00000169174", "ENSG00000141510"]
annotations = gp.convert(organism='hsapiens', query=gene_list)
print(annotations)

gene_set = ["ENSG00000157764", "ENSG00000169174", "ENSG00000141510"]
enrichment_results = gp.profile(organism='hsapiens', query=gene_set)
print(enrichment_results)

# 提取富集分析结果的关键信息
terms = enrichment_results[['native', 'name', 'p_value']]
terms['p_value'] = -1 * terms['p_value'].apply(np.log10)  # 取-p值的对数

# 创建热图
heatmap_data = terms.pivot(index='name', columns='native', values='p_value')
plt.figure(figsize=(10, 8))
sns.heatmap(heatmap_data, cmap='YlGnBu', linewidths=0.5, annot=True, fmt=".2f", cbar=True)
plt.title('Gene Set Enrichment Analysis')
plt.xlabel('Term ID')
plt.ylabel('Term Name')
plt.show()
"""

'''
# compute LR scores and plot interactions among multi-regions
LR_pairs = find_all_LR_pair_expressed(adata, Ligand, Receptor)
region_list = ['IDC_4', 'IDC_2', 'IDC_5', 'IDC_6', 'IDC_8', 'DCIS/LCIS_4', 'DCIS/LCIS_5', 'Tumor_edge_1', 'Tumor_edge_2']
weight = pd.DataFrame(columns=region_list, index=region_list)
for i in range(len(region_list)):
    for j in range(len(region_list)):
        if region_list[i] == region_list[j]:
            weight.loc[region_list[i], region_list[j]] = 0
        else:
            lr_score = compute_LR_score(adata, LR_pairs, lig_cell=region_list[i], rec_cell=region_list[j])
            weight.loc[region_list[i], region_list[j]] = np.sum(lr_score) / lr_score.shape[0]




weight.to_csv('C:/Users/haiyu/Desktop/CCI/weight_new.csv')

for x in region_list:
    adata_x = adata_pearlST[adata_pearlST.obs['annotation'] == x]
    print(x, "The number of cells:", len(adata_x))
'''
