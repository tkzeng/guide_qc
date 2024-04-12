
python guide_assignment_qc.py \
    --guide_h5ad example_data/guide_counts.h5ad \
    --cell_id bc_wells \
    --out_dir example_data/figures \
    --prefix demo \
    --gex_metadata example_data/gex_metadata.tsv \
    --umi_id tscp_count \
    --group_by cluster
