// CPG: DROP results for one (annotation, DROP group) -> seqr RNA-seq upload files; see bin/export_seqr_rna.py
process EXPORT_SEQR {
    tag "${meta.id}"
    label 'process_single'

    container "${params.cpg_container}"

    input:
    tuple val(meta), path(outrider), path(fraser), path(tpm), path(gene_name_mapping)
    path sg_map
    path seqr_gene_ids
    val cohort

    output:
    path 'seqr_export/**', emit: files

    script:
    def inputs = [
        outrider ? "--outrider ${outrider}" : '',
        fraser ? "--fraser ${fraser}" : '',
        tpm ? "--tpm ${tpm}" : '',
        gene_name_mapping ? "--gene-name-mapping ${gene_name_mapping}" : '',
        seqr_gene_ids ? "--seqr-gene-ids ${seqr_gene_ids}" : '',
        params.seqr_project_override ? "--seqr-project-override ${params.seqr_project_override}" : '',
        params.ignore_exome_check ? '--ignore-exome-check' : '',
    ].findAll { arg -> arg }.join(' ')
    """
    export_seqr_rna.py \\
        --sg-map ${sg_map} \\
        --cohort ${cohort} \\
        --ar-guid ${params.arGuid ?: 'local'} \\
        ${inputs} \\
        --output-dir seqr_export

    # one report per (annotation, DROP group), so several tasks can publish into the same folder
    mv seqr_export/export_report.tsv seqr_export/export_report.${meta.id}.tsv
    mv seqr_export/export_skipped.tsv seqr_export/export_skipped.${meta.id}.tsv
    """

    stub:
    """
    mkdir -p seqr_export
    printf 'data_type\\tseqr_project_guid\\tseqr_tissue\\tseqr_product\\tfile\\trows\\tsamples\\n' > seqr_export/export_report.${meta.id}.tsv
    printf 'sg_id\\treason\\n' > seqr_export/export_skipped.${meta.id}.tsv
    """
}