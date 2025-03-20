#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameter validation
def validateParameters() {
    if (!params.reference) {
        exit 1, "Reference genome file not specified!"
    }
    if (!params.fasta_dir) {
        exit 1, "Input FASTA directory not specified!"
    }
}

// Import modules
include { EXTRACT_REFERENCE     } from "${projectDir}/modules/local/extract_reference/main"
include { GOFASTA_VARIANTS      } from "${projectDir}/modules/local/gofasta_variants/main"
include { PULL_DMS              } from "${projectDir}/modules/local/pull_dms/main"
include { INTRAHOST_DMS         } from "${projectDir}/modules/local/intrahost_dms/main"

/*
========================================================================================
    WORKFLOWS
========================================================================================
*/

workflow {
    // Print pipeline info
    log.info """
    ==============================================
    H5N1 DMS ANALYSIS PIPELINE
    ==============================================
    Reference         : ${params.reference}
    Output directory  : ${params.outdir}
    FASTA directory   : ${params.fasta_dir}
    Segments          : ${params.segments.join(', ')}
    ==============================================
    """

    // Validate parameters
    validateParameters()
    
    // Create channel for segments
    Channel
        .fromList(params.segments)
        .set { segments_ch }
    
    // Extract reference segments
    EXTRACT_REFERENCE(segments_ch, params.reference)
    
    // Find consensus sequences by segment
    consensus_by_segment = EXTRACT_REFERENCE.out.reference_segments
        .map { segment, ref_fasta ->
            tuple(
                segment,
                ref_fasta,
                file("${params.fasta_dir}/*_${segment}_cns.fa")
            )
        }
    
    // Find protein variants
    GOFASTA_VARIANTS(consensus_by_segment)

    // Process DMS data
    PULL_DMS(GOFASTA_VARIANTS.out.aa_changes)
    
    // Combine all DMS files
    PULL_DMS.out.dms_data
        .map { it[1] }  // Extract just the JSON files
        .collect()      // Collect all segments
        .set { all_dms_files }
    
    // Process intrahost data with DMS scores
    INTRAHOST_DMS(all_dms_files)
}