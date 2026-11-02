// CPG: Metamist cohort -> DROP sample sheet + SG map; see bin/build_drop_samplesheet.py
process BUILD_DROP_SAMPLESHEET {
    tag "${cohort}"
    label 'process_single'

    container "${params.cpg_container}"

    input:
    val cohort
    val seqr_data_type

    output:
    path 'samplesheet.tsv', emit: samplesheet
    path 'sg_map.tsv'     , emit: sg_map

    script:
    """
    build_drop_samplesheet.py \\
        --cohort ${cohort} \\
        --seqr-data-type ${seqr_data_type} \\
        --output-dir .
    """

    stub:
    """
    printf 'RNA_ID\\tRNA_BAM_FILE\\tRNA_BAI_FILE\\tDNA_VCF_FILE\\tDNA_TBI_FILE\\tDNA_ID\\tDROP_GROUP\\tPAIRED_END\\tCOUNT_MODE\\tCOUNT_OVERLAPS\\tSTRAND\\tHPO_TERMS\\tGENE_COUNTS_FILE\\tGENE_ANNOTATION\\tGENOME\\tSPLICE_COUNTS_DIR\\n' > samplesheet.tsv
    printf 'sg_id\\tparticipant_eid\\tmetamist_project\\tseqr_project_guid\\tsample_type\\tsg_type\\tseqr_tissue\\tseqr_product\\tin_exome_project\\n' > sg_map.tsv
    """
}