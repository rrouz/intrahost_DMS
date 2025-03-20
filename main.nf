#!/usr/bin/env nextflow

include { EXTRACT_REFERENCE } from "${projectDir}/modules/local/extract_reference"
include { GOFASTA_VARIANTS } from "${projectDir}/modules/local/gofasta_variants"
include { PULL_DMS } from "${projectDir}/modules/local/pull_dms"
include { INTRAHOST_DMS } from "${projectDir}/modules/local/intrahost_dms"

workflow {
    segments = Channel.fromList(params.segments)
    
    EXTRACT_REFERENCE(segments, params.reference)
    
    consensus_by_segment = EXTRACT_REFERENCE.out.reference_segments.map { segment, ref_fasta ->
        tuple(
            segment,
            ref_fasta,
            file("${params.fasta_dir}/*_${segment}_cns.fa")
        )
    }
    
    GOFASTA_VARIANTS(consensus_by_segment)

    PULL_DMS(GOFASTA_VARIANTS.out.aa_changes)
    
    all_dms_files = PULL_DMS.out.dms_data.map { it[1] }.collect()
    
    INTRAHOST_DMS(all_dms_files)
}