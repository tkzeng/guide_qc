import argparse
import pathlib
import time

import numpy as np
import pandas as pd
import seaborn as sns

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm


def guides_per_cell(guide_assign_data):
    """
    Distribution for the number of guides assigned to each cell
    """
    n_guides = guide_assign_data.sum(axis=1)
    bins = np.arange(0, n_guides.max() + 1.5) - 0.5
    sns.histplot(n_guides, bins=bins)
    plt.xlabel("Number of guides per cell")
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig(f"{FULL_PREFIX}.guides_per_cell.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def umis_per_guide(guide_assign_data, guide_cts_data):
    """
    Number of guides against number of UMIs per guide
    """
    x = guide_assign_data.sum(axis=1)
    y = guide_cts_data.sum(axis=1)

    # boxplot if few guides per cell; else, hexbin plot
    if np.ptp(x) < 10:
        sns.boxplot(x=x, y=y)
        plt.xlabel("Number of guides")
        plt.ylabel("Number of guide UMIs")
    else:
        plot = sns.jointplot(x=x, y=y, kind="hex", cmap="viridis", norm=LogNorm())
        plot.ax_marg_x.remove()
        plot.ax_marg_y.remove()
        plt.colorbar()

    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("Number of guides")
    plt.ylabel("Number of guide UMIs")
    plt.savefig(f"{FULL_PREFIX}.umis_per_guide.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def top_guides_count_and_proportion(guide_cts_data, k):
    """
    For low MOI experiments.

    Plot 1: Distribution of guide counts for the top k most abundant guides
    Plot 2: Proportion of guide counts for the most abundant guide = (count)/(count+count for other guides)

    :k: number of top guides to plot
    """
    guide_all_cts = np.sort(guide_cts_data, axis=1)[:, ::-1]

    ##### count #####
    counts, color = [], []
    n = len(guide_all_cts)

    for i in range(k):
        counts.append(guide_all_cts[:, i])
        color += [f"#{i + 1} most abundant guide"] * n

    counts_rest = np.sum(guide_all_cts[:, k:], axis=1)
    counts.append(counts_rest)
    color += [f"All other guides"] * n
    counts = np.concatenate(counts)

    sns.histplot(x=counts, hue=color, bins=200, alpha=0.3)
    plt.xlim(0, np.quantile(counts, 0.99))
    plt.yscale("log")
    plt.xlabel("Guide UMIs per cell")
    plt.ylabel("Count (log scale)")
    plt.savefig(f"{FULL_PREFIX}.top_guides_count.png",
                dpi=300, bbox_inches="tight")
    plt.close()

    guide_max = guide_all_cts[:, 0]
    guide_rest = np.sum(guide_all_cts[:, 1:], axis=1)

    ##### proportion #####
    prop = guide_max / (guide_max + guide_rest)
    x = guide_max[~np.isnan(prop)]
    y = prop[~np.isnan(prop)]
    y = y[x < np.quantile(x, 0.99)]
    x = x[x < np.quantile(x, 0.99)]

    plot = sns.jointplot(x=x, y=y, kind="hex", cmap="viridis", norm=LogNorm())
    plot.ax_marg_x.remove()
    plot.ax_marg_y.remove()

    plt.title("Proportion of guide counts for the most abundant guide")
    plt.colorbar()
    plt.xlabel("Number of UMIs")
    plt.ylabel("Proportion")
    plt.savefig(f"{FULL_PREFIX}.top_guide_proportion.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def top_k_guides_count_and_proportion(guide_assign_data, guide_cts_data, k):
    """
    For experiments with MOI > 1.

    For each of the top k most abundant guides per cell:
    Plot 1: Distribution of guide counts for that guide vs the distribution of guide counts for unassigned guides
    Plot 2: Proportion of guide counts for that guide = (count)/(count+count for unassigned guides)

    :k: number of top guides to plot
    """
    guide_assign_cts = guide_assign_data * guide_cts_data
    guide_other_cts = (1 - guide_assign_data) * guide_cts_data

    guide_assign_cts_sorted = np.sort(guide_assign_cts, axis=1)[:, ::-1]
    guide_other_cts_sum = guide_other_cts.sum(axis=1)

    for i in range(k):
        ##### count #####
        guide_k = guide_assign_cts_sorted[:, i]
        guide_current = guide_k[guide_k > 0]
        guide_rest = guide_other_cts_sum[guide_k > 0]

        counts = np.concatenate([guide_current, guide_rest])
        n = len(guide_current)
        color = [f"#{i + 1} most abundant guide"] * n + ["Unassigned guides"] * n
        sns.histplot(x=counts, hue=color, bins=200, alpha=0.3)
        plt.xlabel("Guide UMIs per cell")
        plt.ylabel("Count (log scale)")
        plt.yscale("log")
        plt.xlim(0, np.quantile(counts, 0.99))
        plt.savefig(f"{FULL_PREFIX}.most_abundant_guide_{i + 1}_count.png",
                    dpi=300, bbox_inches="tight")
        plt.close()

        ##### proportion #####
        prop = guide_current / (guide_current + guide_rest)
        x = guide_current[~np.isnan(prop)]
        y = prop[~np.isnan(prop)]
        y = y[x < np.quantile(x, 0.99)]
        x = x[x < np.quantile(x, 0.99)]

        plot = sns.jointplot(x=x, y=y, kind="hex", cmap="viridis", norm=LogNorm())
        plot.ax_marg_x.remove()
        plot.ax_marg_y.remove()

        plt.title(f"Proportion of guide counts for the #{i + 1} most abundant guide")
        plt.colorbar()
        plt.xlabel("Number of UMIs")
        plt.ylabel("Proportion")
        plt.savefig(f"{FULL_PREFIX}.most_abundant_guide_{i + 1}_proportion.png",
                    dpi=300, bbox_inches="tight")
        plt.close()


def compare_to_background(guide_assign_data, guide_cts_data, k):
    """
    Distribution of guide counts for cells assigned to a guide against
    the distribution of counts for all cells with reads mapped to a guide.

    :k: number of guides to plot
    """
    guides = guide_assign_data.columns
    guide_assign_cts = guide_assign_data * guide_cts_data
    guide_other_cts = (1 - guide_assign_data) * guide_cts_data

    for guide in guides[:k]:
        guide_assign_i = guide_assign_cts.loc[:, guide]
        guide_other_i = guide_other_cts.loc[:, guide]

        guide_current = guide_assign_i[guide_assign_i > 0]
        guide_rest = guide_other_i[guide_assign_i == 0]

        counts = np.concatenate([guide_current, guide_rest])
        color = [f"Guide {guide}"] * len(guide_current) + ["Other cells with the guide"] * len(guide_rest)
        sns.histplot(x=counts, hue=color, bins=200, alpha=0.3)
        plt.xlabel("Guide UMIs per cell")
        plt.ylabel("Count (log scale)")
        plt.yscale("log")
        plt.xlim(0, np.quantile(counts, 0.99))
        plt.savefig(f"{FULL_PREFIX}.{guide}_vs_background.png",
                    dpi=300, bbox_inches="tight")
        plt.close()


def guide_vs_gex_umis(gex_metadata, guide_cts_data, guide_cols, transcript_umis):
    """
    Plot 1: Number of transcriptome UMIs against number of guide UMIs
    Plot 2: Number of transcriptome UMIs against number of unique guide UMIs
    """
    merged = gex_metadata.merge(guide_cts_data, left_index=True, right_index=True)

    x = merged[guide_cols].sum(axis=1)
    y = merged[transcript_umis]

    plot = sns.jointplot(x=x, y=y,
                         kind="hex", cmap="viridis", norm=LogNorm())
    plt.xlabel("# of guide UMIs")
    plt.ylabel("# of transcriptome UMIs")

    plot.ax_marg_x.remove()
    plot.ax_marg_y.remove()
    plt.colorbar()

    plt.savefig(f"{FULL_PREFIX}.guide_vs_gex_umis.png",
                dpi=300, bbox_inches="tight")
    plt.close()


def guide_vs_gex_umis_unique(gex_metadata, guide_assign_data, guide_cols, transcript_umis):
    """
    Number of transcriptome UMIs against number of unique guide UMIs
    """
    merged = gex_metadata.merge(guide_assign_data, left_index=True, right_index=True)

    x = merged[guide_cols].sum(axis=1)
    y = merged[transcript_umis]

    if np.ptp(x) < 10:
        sns.boxplot(x=x, y=y)
    else:
        plot = sns.jointplot(x=x, y=y,
                             kind="hex", cmap="viridis", norm=LogNorm())
        plt.xlabel("# of unique guides")
        plt.ylabel("# of transcriptome UMIs")

        plot.ax_marg_x.remove()
        plot.ax_marg_y.remove()
        plt.colorbar()
        plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))

    plt.savefig(f"{FULL_PREFIX}.guide_vs_gex_umis_unique.png",
                dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--guide_metadata", dest="guide_metadata", required=False,
                        help="tsv with cell_id and guide metadata.")
    parser.add_argument("--gex_metadata", dest="gex_metadata", required=False,
                        help="tsv with cell_id and transcriptome metadata.")
    parser.add_argument("--umi_id", dest="umi_id", required=False,
                        help="Column for the number of transcript UMIs in gex_metadata.")
    parser.add_argument("--top_k", dest="top_k", required=False, default=3,
                        help="How many of the most abundant guides per cell to make plots for.")
    parser.add_argument("--group_by", dest="group_by", required=False,
                        help="List of columns in gex_metadata to use for grouping cells "
                             "(for example, clusters or cell cycle state). "
                             "Separate columns by commas. "
                             "All plots will be remade for each of these groupings.")
    parser.add_argument("--out_dir", dest="out_dir", required=False,
                        help="Output directory.")

    required = parser.add_argument_group("required arguments")
    required.add_argument("--guide_assign", dest="guide_assign", required=True,
                          help="tsv with cell_id and guide assignments.")
    required.add_argument("--guide_counts", dest="guide_counts", required=True,
                          help="tsv with cell_id and guide UMI counts.")
    required.add_argument("--cell_id", dest="cell_id", required=True,
                          help="Column for the cell ID.")
    required.add_argument("--prefix", dest="prefix", required=True,
                          help="Prefix for output files.")

    args = parser.parse_args()

    global FULL_PREFIX

    CELL_ID = args.cell_id
    TRANSCRIPT_UMIS = args.umi_id
    TOP_K = args.top_k
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

    guide_assign = pd.read_csv(args.guide_assign, sep='\t')
    guide_cts = pd.read_csv(args.guide_counts, sep='\t')

    GUIDE_COLS = [col for col in guide_assign.columns if col != CELL_ID]
    N_GUIDES = len(GUIDE_COLS)

    guide_assign_data = guide_assign[GUIDE_COLS]
    guide_assign_data.index = guide_assign[CELL_ID]
    guide_cts_data = guide_cts[GUIDE_COLS]
    guide_cts_data.index = guide_cts[CELL_ID]

    # make plots for all cells
    print("making plots for all cells")
    guides_per_cell(guide_assign_data)
    umis_per_guide(guide_assign_data, guide_cts_data)
    top_guides_count_and_proportion(guide_cts_data, k=min(3, N_GUIDES - 1))
    top_k_guides_count_and_proportion(guide_assign_data, guide_cts_data, k=min(TOP_K, N_GUIDES - 1))
    compare_to_background(guide_assign_data, guide_cts_data, k=min(10, N_GUIDES))

    if args.gex_metadata is not None:
        gex_metadata = pd.read_csv(args.gex_metadata, sep='\t')
        gex_metadata.index = gex_metadata[CELL_ID]
        guide_vs_gex_umis(gex_metadata, guide_cts_data, GUIDE_COLS, TRANSCRIPT_UMIS)
        guide_vs_gex_umis(gex_metadata, guide_assign_data, GUIDE_COLS, TRANSCRIPT_UMIS)

        # make plots for groups of cells
        for col in args.group_by.strip().split(','):
            groups = set(gex_metadata[col])
            for group in groups:
                cells = gex_metadata[gex_metadata[col] == group].index
                print(f"making plots for cell group {col}:{group} ({len(cells)} cells)")
                guide_assign_data_sub = guide_assign_data.loc[cells]
                guide_cts_data_sub = guide_cts_data.loc[cells]
                gex_metadata_sub = gex_metadata.loc[cells]

                pathlib.Path(f"{OUT_DIR}/{col}.{group}").mkdir(parents=True, exist_ok=True)
                FULL_PREFIX = f"{OUT_DIR}/{col}.{group}/{PREFIX}"

                guides_per_cell(guide_assign_data_sub)
                umis_per_guide(guide_assign_data_sub, guide_cts_data_sub)
                top_guides_count_and_proportion(guide_cts_data_sub, k=min(3, N_GUIDES - 1))
                top_k_guides_count_and_proportion(guide_assign_data_sub, guide_cts_data_sub, k=TOP_K)
                compare_to_background(guide_assign_data_sub, guide_cts_data_sub, k=min(10, N_GUIDES - 1))
                guide_vs_gex_umis(gex_metadata_sub, guide_cts_data_sub, GUIDE_COLS, TRANSCRIPT_UMIS)
                guide_vs_gex_umis(gex_metadata_sub, guide_assign_data_sub, GUIDE_COLS, TRANSCRIPT_UMIS)
