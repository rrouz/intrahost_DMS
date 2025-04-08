process INTRAHOST_DMS {
    tag "intrahost_analysis"
    label 'process_high'
    publishDir "${params.outdir}/intrahost_dms", mode: params.publish_dir_mode

    input:
    path dms_files

    output:
    path "combined_variants.tsv", emit: combined_variants

    script:
    def args = task.ext.args ?: ''
    """
    intrahost_dms.py \\
        --dms_file '${params.dms_file}' \\
        --data_dir ${params.variants_dir} \\
        --output_dir ./ \\
        --threads ${task.cpus} \\
        --chunk_size 10000 \\
        ${args}
    """
}