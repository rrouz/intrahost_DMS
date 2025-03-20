process GOFASTA_VARIANTS {
    publishDir "${params.outdir}/mutations", mode: 'copy'
    
    input:
    tuple val(segment), path(combined_fasta), path(reference_fasta)
    
    output:
    tuple val(segment), path("aa_changes_${segment}.csv"), emit: aa_changes
    
    script:
    """
    minimap2 -t ${task.cpus} -a -x asm20 --score-N=0 ${reference_fasta} ${combined_fasta} > aligned_${segment}.sam
    gofasta sam variants -s aligned_${segment}.sam  --annotation "${params.data_dir}/gff/${segment}.gff" -o aa_changes_${segment}.csv
    """
}