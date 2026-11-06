#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CPG: build the DROP sample sheet for a Metamist RNA cohort
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Run before main.nf, whose --input is the samplesheet.tsv written here, and whose --sg_map is
    the sg_map.tsv. A separate entry point because DROP parses --input during initialisation.
    See docs/CPG_METAMIST_SEQR_DESIGN.md (T2.3).

    nextflow run prepare.nf --cohort COH0001 --outdir gs://.../drop/COH0001/prepare
----------------------------------------------------------------------------------------
*/

include { BUILD_DROP_SAMPLESHEET } from './modules/local/build_drop_samplesheet'

workflow {
    main:
    if (!params.cohort) {
        error "--cohort is required: the Metamist cohort ID to build the sample sheet for"
    }
    if (!params.outdir) {
        error "--outdir is required: where samplesheet.tsv and sg_map.tsv are published"
    }
    if (!(params.seqr_data_type in ['exome', 'genome'])) {
        error "--seqr_data_type must be exome or genome, got '${params.seqr_data_type}'"
    }

    BUILD_DROP_SAMPLESHEET(params.cohort, params.seqr_data_type)
}