# Guide assignment QC

Makes QC plots for assignments of guides to cells. See `guide_qc.sh` for an example and `guide_assignment_qc.py` for descriptions of all the plots.

Input files need to be generated from previously computed guide assignments. See `parse.scanpy.preprocessing.ipynb` for one way to do this using outputs from the Parse pipeline.

### Installation

Required packages can be installed as follows:
```
pip install numpy pandas seaborn matplotlib 
```

### Usage
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
