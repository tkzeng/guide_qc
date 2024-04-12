import pandas as pd
import scipy
import argparse
import scanpy as sc
import numpy as np


def process_10x(mtx, guide_calls):
    guide_counts = sc.read_10x_mtx(mtx, gex_only=False)
    guide_counts = guide_counts[:, guide_counts.var["feature_types"] == "CRISPR Guide Capture"]

    guide_calls = pd.read_csv(guide_calls, index_col="cell_barcode")
    guide_calls["num_umis"] = guide_calls["num_umis"].apply(
        lambda x: np.sum([int(count) for count in x.split('|')]))
    guide_counts.obs = guide_counts.obs.merge(
        guide_calls, left_index=True, right_index=True, how="outer")

    guide_counts.obs["gene_count"] = guide_counts.obs["num_features"].fillna(0)
    guide_counts.obs["tscp_count"] = guide_counts.obs["num_umis"].fillna(0)
    guide_counts.obs = guide_counts.obs.drop(columns=["num_features", "num_umis"])

    # order guides as in count matrix
    guides = guide_counts.var.index
    guide_to_idx = {guide: i for i, guide in enumerate(guides)}

    guide_calls["guides"] = guide_calls["feature_call"].apply(lambda x: x.split('|'))
    row_idx, col_idx, data = [], [], []

    for i, cell in enumerate(guide_counts.obs.index):
        if cell not in guide_calls.index:
            continue
        for guide in guide_calls.loc[cell]["guides"]:
            row_idx.append(i)
            col_idx.append(guide_to_idx[guide])
            data.append(1)

    guide_assign = scipy.sparse.csr_matrix((data, (row_idx, col_idx)),
                                           dtype=np.int32)
    guide_counts.layers["assignment"] = guide_assign

    return guide_counts


def process_parse(metadata, mtx, guides, guide_calls):
    metadata = pd.read_csv(metadata)
    guide_calls = pd.read_csv(guide_calls)

    cells = metadata["bc_wells"]
    cell_to_idx = {cell: i for i, cell in enumerate(cells)}

    guide_counts = sc.read_mtx(mtx)
    guide_data = pd.read_csv(guides)

    # find genes with nan values and filter
    gene_data = guide_data[guide_data.gene_name.notnull()]
    notNa = gene_data.index
    notNa = notNa.to_list()

    # remove genes with nan values and assign gene names
    guide_counts = guide_counts[:, notNa]
    guide_counts.var = gene_data
    guide_counts.var.set_index("gene_name", inplace=True)
    guide_counts.var.index.name = None
    guide_counts.var_names_make_unique()

    # add cell meta data to anndata object
    guide_counts.obs = metadata
    guide_counts.obs.set_index("bc_wells", inplace=True)
    guide_counts.obs.index.name = None
    guide_counts.obs_names_make_unique()

    # order guides as in count matrix
    guides = guide_counts.var.index
    guide_to_idx = {guide: i for i, guide in enumerate(guides)}

    row_idx, col_idx, data = [], [], []
    for i, row in guide_calls.iterrows():
        row_idx.append(cell_to_idx[row["bc_wells"]])
        col_idx.append(guide_to_idx[row["guide"]])
        data.append(1)

    guide_assign = scipy.sparse.csr_matrix((data, (row_idx, col_idx)),
                                           shape=(len(cells), len(guides)),
                                           dtype=np.int32)
    guide_counts.layers["assignment"] = guide_assign

    return guide_counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    required.add_argument("--guide_calls", dest="guide_calls", required=True,
                          help="10x: tsv formatted as protospacer_calls_per_cell.csv "
                               "(https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-crispr-overview).\n"
                               "Parse: tsv formatted as guide_assignment.csv "
                               "(https://support.parsebiosciences.com/hc/en-us/articles/17166220335636-Pipeline-Setup-and-Use-Current-Version-).")
    required.add_argument("--mtx", dest="mtx", required=True,
                          help="Sparse matrix (Matrix Market format) with guide counts.")
    required.add_argument("--guides", dest="guides", required=False,
                          help="List of guides (required for Parse).")
    required.add_argument("--metadata", dest="metadata", required=True,
                          help="csv with metadata information per cell.")
    required.add_argument("--out_dir", dest="out_dir", required=True,
                          help="Output directory.")
    required.add_argument("--lib_type", choices=["10x", "Parse"], required=True,
                          help="Library type (10x or Parse).")
    args = parser.parse_args()

    if args.lib_type == "10x":
        guide_counts = process_10x(args.mtx, args.guide_calls)
    elif args.lib_type == "Parse":
        guide_counts = process_parse(args.metadata, args.mtx, args.guides, args.guide_calls)

    guide_counts.write_h5ad(f"{args.out_dir}/guide_counts.h5ad")
