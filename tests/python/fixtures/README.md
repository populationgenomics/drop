Outputs of unmodified nf-core/drop `-profile test` (14 public 1000 Genomes samples, chr21, GENCODE v29 lift37),
captured on 2026-11-02 for the seqr export tests (docs/CPG_METAMIST_SEQR_DESIGN.md, T0.4):

- `OUTRIDER_results.tsv`: `processed_results/aberrant_expression/v29/outrider/outrider/`
- `total_counts.Rds`: same folder; input to the TPM step
- `results_per_junction.tsv`: `processed_results/aberrant_splicing/results/v29/fraser/fraser/`
- `gene_name_mapping.tsv`: PREPROCESSGENEANNOTATION's output (not published upstream). Comma-separated despite the name.

The test profile sets `ae_padj_cutoff = 1`, so OUTRIDER_results.tsv holds every gene/sample pair, not just outliers.
