# minimal working script
#python guide_assignment_qc.py \
#    --guide_assign parse/guide_assignment.tsv \
#    --guide_counts parse/guide_counts.tsv \
#    --cell_id bc_wells \
#    --prefix parse_demo

# to make all plots
python guide_assignment_qc.py \
    --guide_assign parse/guide_assignment.tsv \
    --guide_counts parse/guide_counts.tsv \
    --cell_id bc_wells \
    --out_dir qc_figures \
    --prefix parse_demo \
    --guide_metadata parse/guide_metadata.tsv \
    --gex_metadata parse/gex_metadata.tsv \
    --umi_id tscp_count \
    --group_by leiden,pct_counts_mt_below_7

