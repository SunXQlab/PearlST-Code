
import numpy as np
import scanpy as sc
import anndata
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt


def prepare_figure(rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Prepare the figure and axes given the configuration
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    fig, axs = plt.subplots(1, 1, figsize=(csz, rsz))
    plt.subplots_adjust(wspace=wspace, hspace=hspace, left=left, right=right, bottom=bottom, top=top)
    return fig, axs

def prepare_figure_multi_4(rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Prepare the figure and axes given the configuration
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    fig, axs = plt.subplots(1, 4, figsize=(csz, rsz))
    plt.subplots_adjust(wspace=wspace, hspace=hspace, left=left, right=right, bottom=bottom, top=top)
    return fig, axs


def prepare_figure_multi_5(rsz=4., csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Prepare the figure and axes given the configuration
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    fig, axs = plt.subplots(1, 5, figsize=(csz, rsz))
    plt.subplots_adjust(wspace=wspace, hspace=hspace, left=left, right=right, bottom=bottom, top=top)
    return fig, axs


def pseudo_Spatiotemporal_Map(adata_all, emb_name, n_neighbors=20, resolution=1.0):
    """
    Perform pseudo-Spatiotemporal Map for ST data
    :param pSM_values_save_filepath: the default save path for the pSM values
    :type pSM_values_save_filepath: class:`str`, optional, default: "./pSM_values.tsv"
    :param n_neighbors: The size of local neighborhood (in terms of number of neighboring data
    points) used for manifold approximation. See `https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.neighbors.html` for detail
    :type n_neighbors: int, optional, default: 20
    :param resolution: A parameter value controlling the coarseness of the clustering.
    Higher values lead to more clusters. See `https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html` for detail
    :type resolution: float, optional, default: 1.0
    """
    error_message = "No embedding found, please ensure you have run train() method before calculating pseudo-Spatiotemporal Map!"
    max_cell_for_subsampling = 5000
    try:
        print("Performing pseudo-Spatiotemporal Map")
        adata = anndata.AnnData(adata_all.obsm[emb_name])
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X')
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=resolution)
        """
        adata.obs['leiden'] = pd.Categorical(adata_all.obs['pred_label'])
        """
        sc.tl.paga(adata)
        if adata.shape[0] < max_cell_for_subsampling:
            sub_adata_x = adata.X
        else:
            indices = np.arange(adata.shape[0])
            selected_ind = np.random.choice(indices, max_cell_for_subsampling, False)
            sub_adata_x = adata.X[selected_ind, :]
        sum_dists = distance_matrix(sub_adata_x, sub_adata_x).sum(axis=1)
        adata.uns['iroot'] = np.argmax(sum_dists)
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
        pSM_values = adata.obs['dpt_pseudotime'].to_numpy()
        '''
        save_dir = os.path.dirname(pSM_values_save_filepath)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        np.savetxt(pSM_values_save_filepath, pSM_values, fmt='%.5f', header='', footer='', comments='')
        print(
            f"pseudo-Spatiotemporal Map(pSM) calculation complete, pSM values of cells or spots saved at {pSM_values_save_filepath}!")
        '''
        adata_all.obsm['pSM_values'] = pSM_values
    except NameError:
        print(error_message)
    except AttributeError:
        print(error_message)


def plot_pSM_multi_4(adata, adata1, adata2, adata3, adata4, scatter_sz=1., rsz=4.,
             csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Plot the domain segmentation for ST data in spatial
    :param pSM_figure_save_filepath: the default save path for the figure
    :type pSM_figure_save_filepath: class:`str`, optional, default: "./Spatiotemporal-Map.pdf"
    :param colormap: The colormap to use. See `https://www.fabiocrameri.ch/colourmaps-userguide/` for name list of colormaps
    :type colormap: str, optional, default: roma
    :param scatter_sz: The marker size in points**2
    :type scatter_sz: float, optional, default: 1.0
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    error_message = "No pseudo Spatiotemporal Map data found, please ensure you have run the pseudo_Spatiotemporal_Map() method."
    # try:
    fig, ax = prepare_figure_multi_4(rsz=rsz, csz=csz, wspace=wspace, hspace=hspace, left=left, right=right,
                                      bottom=bottom, top=top)
    x0, y0 = adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]
    ax[0].scatter(x0, y0, s=scatter_sz, c=adata.obsm['pSM_values'], cmap='summer', marker=".")
    ax[0].invert_yaxis()
    # clb = fig.colorbar(st0)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[0].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    ax[0].spines['right'].set_visible(False)  # 去掉边框
    ax[0].spines['top'].set_visible(False)  # 去掉边框
    ax[0].spines['left'].set_visible(False)  # 去掉边框
    ax[0].spines['bottom'].set_visible(False)  # 去掉边框
    ax[0].get_yaxis().set_visible(False)
    ax[0].get_xaxis().set_visible(False)
    # ax[0].set_title("SpaceFlow pSM", fontsize=12)
    ax[0].set_facecolor("none")

    x1, y1 = adata1.obsm["spatial"][:, 0], adata1.obsm["spatial"][:, 1]
    ax[1].scatter(x1, y1, s=scatter_sz, c=adata1.obsm['pSM_values'], cmap='summer', marker=".")
    ax[1].invert_yaxis()
    ax[1].spines['right'].set_visible(False)  # 去掉边框
    ax[1].spines['top'].set_visible(False)  # 去掉边框
    ax[1].spines['left'].set_visible(False)  # 去掉边框
    ax[1].spines['bottom'].set_visible(False)  # 去掉边框
    ax[1].get_yaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st1)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[1].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[1].set_title("STAGATE pSM", fontsize=12)
    ax[1].set_facecolor("none")

    x2, y2 = adata2.obsm["spatial"][:, 0], adata2.obsm["spatial"][:, 1]
    ax[2].scatter(x2, y2, s=scatter_sz, c=adata2.obsm['pSM_values'], cmap='summer', marker=".")
    ax[2].invert_yaxis()
    ax[2].spines['right'].set_visible(False)  # 去掉边框
    ax[2].spines['top'].set_visible(False)  # 去掉边框
    ax[2].spines['left'].set_visible(False)  # 去掉边框
    ax[2].spines['bottom'].set_visible(False)  # 去掉边框
    ax[2].get_yaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st2)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_title("SEDR pSM", fontsize=12)
    ax[2].set_facecolor("none")

    x3, y3 = adata3.obsm["spatial"][:, 0], adata3.obsm["spatial"][:, 1]
    ax[3].scatter(x3, y3, s=scatter_sz, c=adata3.obsm['pSM_values'], cmap='summer', marker=".")
    ax[3].invert_yaxis()
    ax[3].spines['right'].set_visible(False)  # 去掉边框
    ax[3].spines['top'].set_visible(False)  # 去掉边框
    ax[3].spines['left'].set_visible(False)  # 去掉边框
    ax[3].spines['bottom'].set_visible(False)  # 去掉边框
    ax[3].get_yaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st3)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[3].set_title("Scanpy pSM", fontsize=12)
    ax[3].set_facecolor("none")






def plot_pSM_multi_5(adata, adata1, adata2, adata3, adata4, scatter_sz=1., rsz=4.,
             csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Plot the domain segmentation for ST data in spatial
    :param pSM_figure_save_filepath: the default save path for the figure
    :type pSM_figure_save_filepath: class:`str`, optional, default: "./Spatiotemporal-Map.pdf"
    :param colormap: The colormap to use. See `https://www.fabiocrameri.ch/colourmaps-userguide/` for name list of colormaps
    :type colormap: str, optional, default: roma
    :param scatter_sz: The marker size in points**2
    :type scatter_sz: float, optional, default: 1.0
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    error_message = "No pseudo Spatiotemporal Map data found, please ensure you have run the pseudo_Spatiotemporal_Map() method."
    # try:
    fig, ax = prepare_figure_multi_4(rsz=rsz, csz=csz, wspace=wspace, hspace=hspace, left=left, right=right,
                                      bottom=bottom, top=top)
    x0, y0 = adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]
    ax[0].scatter(x0, y0, s=scatter_sz, c=adata.obsm['pSM_values'], cmap='summer', marker=".")
    ax[0].invert_yaxis()
    # clb = fig.colorbar(st0)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[0].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    ax[0].spines['right'].set_visible(False)  # 去掉边框
    ax[0].spines['top'].set_visible(False)  # 去掉边框
    ax[0].spines['left'].set_visible(False)  # 去掉边框
    ax[0].spines['bottom'].set_visible(False)  # 去掉边框
    ax[0].get_yaxis().set_visible(False)
    ax[0].get_xaxis().set_visible(False)
    # ax[0].set_title("SpaceFlow pSM", fontsize=12)
    ax[0].set_facecolor("none")

    x1, y1 = adata1.obsm["spatial"][:, 0], adata1.obsm["spatial"][:, 1]
    ax[1].scatter(x1, y1, s=scatter_sz, c=adata1.obsm['pSM_values'], cmap='summer', marker=".")
    ax[1].invert_yaxis()
    ax[1].spines['right'].set_visible(False)  # 去掉边框
    ax[1].spines['top'].set_visible(False)  # 去掉边框
    ax[1].spines['left'].set_visible(False)  # 去掉边框
    ax[1].spines['bottom'].set_visible(False)  # 去掉边框
    ax[1].get_yaxis().set_visible(False)
    ax[1].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st1)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[1].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[1].set_title("STAGATE pSM", fontsize=12)
    ax[1].set_facecolor("none")

    x2, y2 = adata2.obsm["spatial"][:, 0], adata2.obsm["spatial"][:, 1]
    ax[2].scatter(x2, y2, s=scatter_sz, c=adata2.obsm['pSM_values'], cmap='summer', marker=".")
    ax[2].invert_yaxis()
    ax[2].spines['right'].set_visible(False)  # 去掉边框
    ax[2].spines['top'].set_visible(False)  # 去掉边框
    ax[2].spines['left'].set_visible(False)  # 去掉边框
    ax[2].spines['bottom'].set_visible(False)  # 去掉边框
    ax[2].get_yaxis().set_visible(False)
    ax[2].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st2)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_title("SEDR pSM", fontsize=12)
    ax[2].set_facecolor("none")

    x3, y3 = adata3.obsm["spatial"][:, 0], adata3.obsm["spatial"][:, 1]
    ax[3].scatter(x3, y3, s=scatter_sz, c=adata3.obsm['pSM_values'], cmap='summer', marker=".")
    ax[3].invert_yaxis()
    ax[3].spines['right'].set_visible(False)  # 去掉边框
    ax[3].spines['top'].set_visible(False)  # 去掉边框
    ax[3].spines['left'].set_visible(False)  # 去掉边框
    ax[3].spines['bottom'].set_visible(False)  # 去掉边框
    ax[3].get_yaxis().set_visible(False)
    ax[3].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st3)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[3].set_title("Scanpy pSM", fontsize=12)
    ax[3].set_facecolor("none")

    x4, y4 = adata4.obsm["spatial"][:, 0], adata4.obsm["spatial"][:, 1]
    ax[4].scatter(x4, y4, s=scatter_sz, c=adata4.obsm['pSM_values'], cmap='summer', marker=".")
    ax[4].invert_yaxis()
    ax[4].spines['right'].set_visible(False)  # 去掉边框
    ax[4].spines['top'].set_visible(False)  # 去掉边框
    ax[4].spines['left'].set_visible(False)  # 去掉边框
    ax[4].spines['bottom'].set_visible(False)  # 去掉边框
    ax[4].get_yaxis().set_visible(False)
    ax[4].get_xaxis().set_visible(False)
    # clb = fig.colorbar(st3)
    # clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[2].set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
    # ax[4].set_title("PearlST pSM", fontsize=12)
    ax[4].set_facecolor("none")



    '''
        save_dir = os.path.dirname(pSM_figure_save_filepath)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(pSM_figure_save_filepath, dpi=300)
        print(f"Plotting complete, pseudo-Spatiotemporal Map figure saved at {pSM_figure_save_filepath} !")
        plt.close('all')
    '''
    # except NameError:
    #     print(error_message)
    # except AttributeError:
    #     print(error_message)


def plot_pSM(adata, scatter_sz=1., rsz=4.,
             csz=4., wspace=.4, hspace=.5, left=0.125, right=0.9, bottom=0.1, top=0.9):
    """
    Plot the domain segmentation for ST data in spatial
    :param pSM_figure_save_filepath: the default save path for the figure
    :type pSM_figure_save_filepath: class:`str`, optional, default: "./Spatiotemporal-Map.pdf"
    :param colormap: The colormap to use. See `https://www.fabiocrameri.ch/colourmaps-userguide/` for name list of colormaps
    :type colormap: str, optional, default: roma
    :param scatter_sz: The marker size in points**2
    :type scatter_sz: float, optional, default: 1.0
    :param rsz: row size of the figure in inches, default: 4.0
    :type rsz: float, optional
    :param csz: column size of the figure in inches, default: 4.0
    :type csz: float, optional
    :param wspace: the amount of width reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type wspace: float, optional
    :param hspace: the amount of height reserved for space between subplots, expressed as a fraction of the average axis width, default: 0.4
    :type hspace: float, optional
    :param left: the leftmost position of the subplots of the figure in fraction, default: 0.125
    :type left: float, optional
    :param right: the rightmost position of the subplots of the figure in fraction, default: 0.9
    :type right: float, optional
    :param bottom: the bottom position of the subplots of the figure in fraction, default: 0.1
    :type bottom: float, optional
    :param top: the top position of the subplots of the figure in fraction, default: 0.9
    :type top: float, optional
    """
    error_message = "No pseudo Spatiotemporal Map data found, please ensure you have run the pseudo_Spatiotemporal_Map() method."
    try:
        fig, ax = prepare_figure(rsz=rsz, csz=csz, wspace=wspace, hspace=hspace, left=left, right=right,
                                      bottom=bottom, top=top)
        x1, y1 = adata.obsm["spatial"][:, 0], adata.obsm["spatial"][:, 1]
        st = ax.scatter(x1, y1, s=scatter_sz, c=adata.obsm['pSM_values'], cmap='summer', marker=".")
        ax.invert_yaxis()
        clb = fig.colorbar(st, shrink=0.4)
        clb.ax.set_ylabel("pseudotime", labelpad=10, rotation=270, fontsize=10, weight='bold')
        ax.set_title("SpaceFlow PSM", fontsize=14)
        ax.set_facecolor("none")
        '''
        save_dir = os.path.dirname(pSM_figure_save_filepath)
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(pSM_figure_save_filepath, dpi=300)
        print(f"Plotting complete, pseudo-Spatiotemporal Map figure saved at {pSM_figure_save_filepath} !")
        plt.close('all')
        '''
    except NameError:
        print(error_message)
    except AttributeError:
        print(error_message)

