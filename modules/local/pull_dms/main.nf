process PULL_DMS {
    tag "${segment}"
    label 'process_medium'
    publishDir "${params.outdir}/dms", mode: params.publish_dir_mode
    
    input:
    tuple val(segment), path(aa_changes_csv)
    
    output:
    tuple val(segment), path("dms_${segment}.tsv"), emit: dms_data
    
    script:
    def args = task.ext.args ?: ''
    """
    pull_dms.py \\
        --dms-file '${params.dms_file}' \\
        --dms-file-2 '${params.dms_file_2}' \\
        --h5_site_header "${params.h5_site_header}" \\
        --mutation_file ${aa_changes_csv} \\
        --output_file dms_${segment}.tsv \\
        ${args}
    """
}