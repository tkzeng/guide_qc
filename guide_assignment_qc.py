import argparse
import pathlib

import numpy as np
import pandas as pd
import seaborn as sns
import scipy

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm

import scanpy as sc


def reads_per_umi(guide_h5ad):
    if "mread_count" not in guide_h5ad.obs:
        return
    x = guide_h5ad.obs["tscp_count"].to_numpy()
    y = guide_h5ad.obs["mread_count"].to_numpy()
    y = y[(x < np.percentile(x, 95)) & (x > 0)]
    x = x[(x < np.percentile(x, 95)) & (x > 0)]
    plot = sns.jointplot(x=x, y=y / x, kind="hex", cmap="viridis", norm=LogNorm(), xscale="log", yscale="log")
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Number of cells", rotation=270)
    plot.ax_marg_x.remove()
    plot.ax_marg_y.remove()
    plt.xlabel("# of UMIs")
    plt.ylabel("# of reads")
    plt.savefig(f"{FULL_PREFIX}.reads_per_umi.png", dpi=300, bbox_inches="tight")
    plt.close()


def cells_per_guide(guide_h5ad):
    cells_per_guide = guide_h5ad.layers["assignment"].sum(axis=0)
    cells_per_guide = np.array(cells_per_guide).flatten()
    sns.histplot(cells_per_guide)
    plt.xlabel("Number of cells per guide")
    plt.savefig(f"{FULL_PREFIX}.cells_per_guide.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def guides_per_cell(guide_h5ad):
    """
    Distribution for the number of guides assigned to each cell
    """
    n_guides = guide_h5ad.layers["assignment"].sum(axis=1)
    n_guides = np.array(n_guides).flatten()
    n_guides = n_guides[n_guides <= np.percentile(n_guides, 95)]
    bins = np.arange(0, n_guides.max() + 1.5) - 0.5
    sns.histplot(n_guides, bins=bins, stat="density")
    plt.xlabel("Number of guides per cell")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

    n75 = int(np.percentile(n_guides, 75))
    n50 = int(np.percentile(n_guides, 50))
    n25 = int(np.percentile(n_guides, 25))
    plt.text(n75 + 0.2, 0.01, "75th percentile", rotation=90)
    plt.text(n50 + 0.2, 0.01, "50th percentile", rotation=90)
    plt.text(n25 + 0.2, 0.01, "25th percentile", rotation=90)
    plt.axvline(n75, linestyle="--", color="black")
    plt.axvline(n50, linestyle="--", color="black")
    plt.axvline(n25, linestyle="--", color="black")

    x, y, hue = [], [], []
    for lambda_param in [1, 3, 5, 8, 11]:
        k_values = np.arange(0, np.max(n_guides))
        pmf_values = scipy.stats.poisson.pmf(k_values, lambda_param)
        x += k_values.tolist()
        y += pmf_values.tolist()
        hue += ([lambda_param] * len(k_values))
    sns.lineplot(x=x, y=y, hue=hue)
    plt.legend(title="Poisson λ")

    plt.xlim(0, int(np.percentile(n_guides, 95)))
    plt.savefig(f"{FULL_PREFIX}.guides_per_cell.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    return n75


def umis(guide_h5ad):
    """
    Plot 1: Number of UMIs per cell (histogram)
    Plot 2: Number of UMIs per guide (histogram)
    """
    x = guide_h5ad.layers["assignment"].sum(axis=1)
    x = np.array(x).flatten()
    y = guide_h5ad.X.sum(axis=1)
    y = np.array(y).flatten()

    ### UMIs per cell ###

    n75 = int(np.percentile(y, 75))
    n50 = int(np.percentile(y, 50))
    n25 = int(np.percentile(y, 25))
    plt.text(n75 + 0.2, 0.01, "75th percentile", rotation=90)
    plt.text(n50 + 0.2, 0.01, "50th percentile", rotation=90)
    plt.text(n25 + 0.2, 0.01, "25th percentile", rotation=90)
    plt.axvline(n75, linestyle="--", color="black")
    plt.axvline(n50, linestyle="--", color="black")
    plt.axvline(n25, linestyle="--", color="black")

    sns.histplot(y, log_scale=False)
    plt.xlim(0, np.percentile(y, 95))
    plt.xlabel("Number of UMIs per cell")
    plt.savefig(f"{FULL_PREFIX}.umis_per_cell.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    ### UMIs per guide ###

    umis_per_guide = y / x
    umis_per_guide[~np.isfinite(umis_per_guide)] = 0

    n75 = int(np.percentile(umis_per_guide, 75))
    n50 = int(np.percentile(umis_per_guide, 50))
    n25 = int(np.percentile(umis_per_guide, 25))
    plt.text(n75 + 0.2, 0.01, "75th percentile", rotation=90)
    plt.text(n50 + 0.2, 0.01, "50th percentile", rotation=90)
    plt.text(n25 + 0.2, 0.01, "25th percentile", rotation=90)
    plt.axvline(n75, linestyle="--", color="black")
    plt.axvline(n50, linestyle="--", color="black")
    plt.axvline(n25, linestyle="--", color="black")

    sns.histplot(umis_per_guide, log_scale=False)
    plt.xlim(0, np.percentile(umis_per_guide, 95))
    plt.xlabel("Number of UMIs per guide")
    plt.savefig(f"{FULL_PREFIX}.umis_per_guide.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def umis_by_nguides(guide_h5ad):
    """
    Plot 1: Number of guides vs. number of UMIs per cell
    Plot 2: Number of guides vs. number of UMIs per guide
    """
    x = guide_h5ad.layers["assignment"].sum(axis=1)
    x = np.array(x).flatten()
    y = guide_h5ad.X.sum(axis=1)
    y = np.array(y).flatten()

    ### Number of UMIs per cell ###

    box = pd.DataFrame({"Number of guides": x,
                        "Number of UMIs": y})
    plot = sns.boxplot(box,
                       x="Number of guides",
                       y="Number of UMIs", fliersize=0)
    xlim = np.percentile(x, 95)
    ylim = np.percentile(y, 95)

    medians = box.groupby("Number of guides").median().to_numpy().flatten()
    for xtick in plot.get_xticks():
        if xtick < xlim and medians[xtick] < ylim:
            plt.text(xtick, medians[xtick], int(medians[xtick]),
                     horizontalalignment="center", size="x-small",
                     color="black", weight="normal", rotation="vertical")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlim(0, xlim)
    plt.ylim(0, ylim)
    plt.savefig(f"{FULL_PREFIX}.umis_per_cell.n_guides.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    ### Number of UMIs per guide ###

    y1 = y / x
    y1[~np.isfinite(y1)] = 0
    box = pd.DataFrame({"Number of guides": x,
                        "Number of UMIs per guide": y1})
    plot = sns.boxplot(box,
                       x="Number of guides",
                       y="Number of UMIs per guide", fliersize=0)
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlim(0, xlim)
    plt.ylim(0, np.percentile(y1, 95))
    plt.savefig(f"{FULL_PREFIX}.umis_per_guide.n_guides.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def sort_csr_rows_by_value(matrix):
    """
    Sorts each row of a CSR matrix by its values.

    Parameters:
    - matrix: A scipy.sparse.csr_matrix object

    Returns:
    - A new csr_matrix object with each row sorted by value.
    """
    rows, cols = matrix.shape
    sorted_data = []
    sorted_indices = []
    sorted_indptr = [0]

    for i in range(rows):
        row_start = matrix.indptr[i]
        row_end = matrix.indptr[i + 1]
        row_data = matrix.data[row_start:row_end]

        # Sort the row data and get the sorted indices
        sorted_row_indices = np.arange(len(row_data))
        sorted_row_data = np.sort(row_data)[::-1]

        # Append the sorted row to the new data and indices lists
        sorted_data.extend(sorted_row_data)
        sorted_indices.extend(sorted_row_indices)
        sorted_indptr.append(sorted_indptr[-1] + len(row_data))

    # Create a new CSR matrix with the sorted data
    sorted_matrix = scipy.sparse.csr_matrix((sorted_data, sorted_indices, sorted_indptr), shape=(rows, cols))
    return sorted_matrix


def top_k_guides_count_and_proportion(guide_h5ad, k):
    """
    For each of the top k most abundant assigned guides per cell:
    Plot 1: Distribution of guide counts for that guide vs the distribution of guide counts for ALL unassigned guides
    Plot 2: Proportion of guide counts for that guide = (count)/(count+count for ALL unassigned guides)

    :k: number of top guides to plot
    """
    guide_assign_data = guide_h5ad.layers["assignment"]
    guide_cts_data = guide_h5ad.X

    guide_assign_cts = guide_assign_data.multiply(guide_cts_data)
    guide_assign_cts_sorted = sort_csr_rows_by_value(guide_assign_cts)

    guide_other_cts = guide_cts_data.copy()
    guide_other_cts[guide_assign_data.astype(bool)] = 0
    # guide_other_cts_sorted = sort_csr_rows_by_value(guide_other_cts)

    guide_other_cts_sum = guide_cts_data.sum(axis=1) - guide_assign_cts.sum(axis=1)
    guide_other_cts_sum = np.array(guide_other_cts_sum).flatten()

    for i in range(k):
        ### count ###
        guide_k = guide_assign_cts_sorted[:, i].todense()
        guide_k = np.array(guide_k).flatten()
        guide_current = guide_k[guide_k > 0]
        guide_rest = guide_other_cts_sum[guide_k > 0]
        # guide_rest = guide_other_cts_sorted[:,0].todense()
        # guide_rest = np.array(guide_rest).flatten()[guide_k > 0]

        counts = np.concatenate([guide_current, guide_rest])
        n = len(guide_current)
        color = [f"#{i + 1} most abundant guide"] * n + ["Unassigned guides"] * n
        color = np.array(color)[counts <= np.quantile(guide_current, 0.95)]
        counts = counts[counts <= np.quantile(guide_current, 0.95)]

        n50 = int(np.percentile(guide_current, 50))
        plt.text(n50 + 0.2, 0.01, "50th percentile", rotation=90)
        plt.axvline(n50, linestyle="--", color="black")

        sns.histplot(x=counts, hue=color, alpha=0.3, stat="percent", common_norm=False, log_scale=False)
        plt.xlabel("Guide UMIs per cell")
        plt.ylabel("Percent")
        plt.savefig(f"{FULL_PREFIX}.most_abundant_guide_{i + 1}_count.png",
                    dpi=300, bbox_inches="tight")
        plt.close()

        ### proportion ###
        prop = guide_current / (guide_current + guide_rest)
        x = guide_current[~np.isnan(prop)]
        y = prop[~np.isnan(prop)]
        y = y[x < np.quantile(x, 0.95)]
        x = x[x < np.quantile(x, 0.95)]

        if len(x) == 0:
            continue

        plot = sns.jointplot(x=x, y=y, kind="hex", cmap="viridis", norm=LogNorm())
        plot.ax_marg_x.remove()
        plot.ax_marg_y.remove()

        plt.title(f"Proportion of guide counts for the #{i + 1} most abundant guide")
        cbar = plt.colorbar()
        cbar.ax.set_ylabel("Number of cells", rotation=270)
        plt.xlabel("Number of UMIs")
        plt.ylabel("Proportion = UMIs/(UMIs+UMIs for unassigned guides)")
        plt.savefig(f"{FULL_PREFIX}.most_abundant_guide_{i + 1}_proportion.png",
                    dpi=300, bbox_inches="tight")
        plt.close()


def compare_to_background(guide_h5ad, k):
    """
    Distribution of guide counts for cells assigned to a guide against
    the distribution of counts for all other cells with reads mapped to a guide.

    :k: number of guides to plot
    """
    guide_assign_data = guide_h5ad.layers["assignment"]
    guide_cts_data = guide_h5ad.X

    guide_assign_cts = guide_assign_data.multiply(guide_cts_data)
    guide_other_cts = guide_cts_data.copy()
    guide_other_cts[guide_assign_data.astype(bool)] = 0

    guides = guide_h5ad.var.index

    for i in range(k):
        guide_assign_i = guide_assign_cts[:, i].todense()
        guide_other_i = guide_other_cts[:, i].todense()

        guide_current = guide_assign_i[guide_assign_i > 0]
        guide_current = np.array(guide_current).flatten()
        if len(guide_current) == 0:
            continue
        guide_rest = guide_other_i[guide_assign_i == 0]
        guide_rest = np.array(guide_rest).flatten()

        counts = np.concatenate([guide_current, guide_rest])
        color = [f"{guides[i]}"] * len(guide_current) + ["Other cells with the guide"] * len(guide_rest)
        color = np.array(color)[counts <= np.quantile(guide_current, 0.95)]
        counts = counts[counts <= np.quantile(guide_current, 0.95)]

        n50 = int(np.percentile(guide_current, 50))
        plt.text(n50 + 0.2, 0.01, "50th percentile", rotation=90)
        plt.axvline(n50, linestyle="--", color="black")

        sns.histplot(x=counts, hue=color, alpha=0.3, log_scale=False, stat="percent", common_norm=False)
        plt.xlabel("Guide UMIs per cell")
        plt.ylabel("Percent")
        plt.savefig(f"{FULL_PREFIX}.{guides[i]}_vs_background.png",
                    dpi=300, bbox_inches="tight")
        plt.close()


def guide_vs_gex_umis(gex_metadata, guide_h5ad, tscp_col):
    """
    Plot 1: Number of transcriptome UMIs against number of guide UMIs
    Plot 2: Number of transcriptome UMIs against number of unique guide UMIs
    """
    idx = list(set(guide_h5ad.obs.index) & set(gex_metadata.index))
    x = guide_h5ad.obs.loc[idx, "tscp_count"]
    y = gex_metadata.loc[idx, tscp_col]
    y = y[x < np.quantile(x, 0.95)]
    x = x[x < np.quantile(x, 0.95)]

    plot = sns.jointplot(x=x, y=y,
                         kind="hex", cmap="viridis", norm=LogNorm())
    plt.xlabel("# of guide UMIs")
    plt.ylabel("# of transcriptome UMIs")
    plot.ax_marg_x.remove()
    plot.ax_marg_y.remove()
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Number of cells", rotation=270)
    plt.savefig(f"{FULL_PREFIX}.guide_vs_gex_umis.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    x = guide_h5ad.obs.loc[idx, "gene_count"]
    y = gex_metadata.loc[idx, tscp_col]
    y = y[x < np.quantile(x, 0.95)]
    x = x[x < np.quantile(x, 0.95)]

    plot = sns.jointplot(x=x, y=y,
                         kind="hex", cmap="viridis", norm=LogNorm())
    plt.xlabel("# of guide UMIs")
    plt.ylabel("# of transcriptome UMIs")
    plot.ax_marg_x.remove()
    plot.ax_marg_y.remove()
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Number of cells", rotation=270)
    plt.savefig(f"{FULL_PREFIX}.guide_vs_gex_umis_unique.png",
                dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--gex_metadata", dest="gex_metadata", required=False,
                        help="tsv with cell_id and transcriptome metadata.")
    parser.add_argument("--umi_id", dest="umi_id", required=False,
                        help="Column for the number of transcript UMIs in gex_metadata.")
    parser.add_argument("--group_by", dest="group_by", required=False,
                        help="List of columns in gex_metadata to use for grouping cells "
                             "(for example, clusters or cell cycle state). "
                             "Separate columns by commas. "
                             "All plots will be remade for each of these groupings.")
    parser.add_argument("--out_dir", dest="out_dir", required=False,
                        help="Output directory.")

    required = parser.add_argument_group("required arguments")
    parser.add_argument("--guide_h5ad", dest="guide_h5ad", required=True)
    required.add_argument("--cell_id", dest="cell_id", required=True,
                          help="Column for the cell ID.")
    required.add_argument("--prefix", dest="prefix", required=True,
                          help="Prefix for output files.")

    args = parser.parse_args()

    plt.rcParams["figure.figsize"] = (5, 5)

    global FULL_PREFIX

    CELL_ID = args.cell_id
    TRANSCRIPT_UMIS = args.umi_id
    PREFIX = args.prefix
    if args.out_dir is None:
        OUT_DIR = './'
    else:
        OUT_DIR = args.out_dir
    FULL_PREFIX = f"{OUT_DIR}/{PREFIX}"

    # make output directory
    pathlib.Path(OUT_DIR).mkdir(parents=True, exist_ok=True)

    flog = open(f"{FULL_PREFIX}.log", 'w')
    print(vars(args), file=flog)

    guide_h5ad = sc.read_h5ad(args.guide_h5ad)
    N_GUIDES = len(guide_h5ad.var)

    print("making plots for all cells")

    n75 = guides_per_cell(guide_h5ad)
    cells_per_guide(guide_h5ad)
    reads_per_umi(guide_h5ad)
    umis(guide_h5ad)
    umis_by_nguides(guide_h5ad)
    top_k_guides_count_and_proportion(guide_h5ad, k=n75)
    compare_to_background(guide_h5ad, k=min(25, N_GUIDES))

    if args.gex_metadata is not None:
        gex_metadata = pd.read_csv(args.gex_metadata, sep='\t')
        gex_metadata.index = gex_metadata[CELL_ID]
        guide_vs_gex_umis(gex_metadata, guide_h5ad, args.umi_id)

        if args.group_by is not None:
            for col in args.group_by.strip().split(','):
                groups = set(gex_metadata[col])
                for group in groups:
                    cells = gex_metadata[gex_metadata[col] == group].index
                    print(f"making plots for cell group {col}:{group} ({len(cells)} cells)")
                    guide_h5ad_sub = guide_h5ad[cells]
                    gex_metadata_sub = gex_metadata.loc[cells]

                    pathlib.Path(f"{OUT_DIR}/{col}.{group}").mkdir(parents=True, exist_ok=True)
                    FULL_PREFIX = f"{OUT_DIR}/{col}.{group}/{PREFIX}"

                    n75 = guides_per_cell(guide_h5ad_sub)
                    cells_per_guide(guide_h5ad)
                    reads_per_umi(guide_h5ad_sub)
                    umis(guide_h5ad_sub)
                    umis_by_nguides(guide_h5ad_sub)
                    top_k_guides_count_and_proportion(guide_h5ad_sub, k=n75)
                    compare_to_background(guide_h5ad_sub, k=min(25, N_GUIDES))
                    guide_vs_gex_umis(gex_metadata, guide_h5ad_sub, args.umi_id)
