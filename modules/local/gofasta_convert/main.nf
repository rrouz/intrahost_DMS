process GOFASTA_CONVERT {
    publishDir "${params.outdir}/mutations", mode: 'copy', pattern: '*.json'
    
    input:
    tuple val(segment), path(aa_changes_csv)
    
    output:
    tuple val(segment), path("aa_changes_${segment}.json"), emit: aa_changes_json
    tuple val(segment), path("aa_changes_${segment}.csv"), emit: aa_changes_csv
    
    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/bin/convert_mutations.py \\
        --input "${aa_changes_csv}" \\
        --output "aa_changes_${segment}.json" \\
        ${args}
    """
} 