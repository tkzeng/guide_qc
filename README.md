# Guide assignment QC

Makes QC plots for assignments of guides to cells. See `guide_qc.sh` for an example and `guide_assignment_qc.py` for descriptions of all the plots.

Input files need to be generated from previously computed guide assignments. See `preprocess.py` for one way to do this using outputs from the 10x or Parse pipelines.

### Installation

Required packages can be installed as follows:
```
pip install numpy pandas seaborn matplotlib scipy 'scanpy[leiden]'
```

### Usage
Script to preprocess the outputs of the standard 10x or Parse pipelines:
```
usage: preprocess.py [-h] --guide_calls GUIDE_CALLS --mtx MTX [--guides GUIDES] --metadata METADATA --out_dir OUT_DIR --lib_type {10x,Parse}

options:
  -h, --help            show this help message and exit

required arguments:
  --guide_calls GUIDE_CALLS
                        10x: tsv formatted as protospacer_calls_per_cell.csv (https://www.10xgenomics.com/support/software/cell-
                        ranger/latest/analysis/outputs/cr-outputs-crispr-overview). Parse: tsv formatted as guide_assignment.csv
                        (https://support.parsebiosciences.com/hc/en-us/articles/17166220335636-Pipeline-Setup-and-Use-Current-Version-). (default: None)
  --mtx MTX             Sparse matrix (Matrix Market format) with guide counts. (default: None)
  --guides GUIDES       List of guides (required for Parse). (default: None)
  --metadata METADATA   csv with metadata information per cell. (default: None)
  --out_dir OUT_DIR     Output directory. (default: None)
  --lib_type {10x,Parse}
                        Library type (10x or Parse). (default: None)
```

Script to QC guide assignment:
```
usage: guide_assignment_qc.py [-h] [--guide_metadata GUIDE_METADATA] [--gex_metadata GEX_METADATA] [--umi_id UMI_ID] [--top_k TOP_K]
                              [--group_by GROUP_BY] [--out_dir OUT_DIR] --guide_assign GUIDE_ASSIGN --guide_counts GUIDE_COUNTS --cell_id CELL_ID
                              --prefix PREFIX

optional arguments:
  -h, --help            show this help message and exit
  --guide_metadata GUIDE_METADATA
                        tsv with cell_id and guide metadata. (default: None)
  --gex_metadata GEX_METADATA
                        tsv with cell_id and transcriptome metadata. (default: None)
  --umi_id UMI_ID       Column for the number of transcript UMIs in gex_metadata. (default: None)
  --top_k TOP_K         How many of the most abundant guides per cell to make plots for. (default: 3)
  --group_by GROUP_BY   List of columns in gex_metadata to use for grouping cells (for example, clusters or cell cycle state). Separate columns by
                        commas. All plots will be remade for each of these groupings. (default: None)
  --out_dir OUT_DIR     Output directory. (default: None)

required arguments:
  --guide_assign GUIDE_ASSIGN
                        tsv with cell_id and guide assignments. (default: None)
  --guide_counts GUIDE_COUNTS
                        tsv with cell_id and guide UMI counts. (default: None)
  --cell_id CELL_ID     Column for the cell ID. (default: None)
  --prefix PREFIX       Prefix for output files. (default: None)
```
