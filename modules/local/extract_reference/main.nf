process EXTRACT_REFERENCE {
    input:
    val segment
    path reference
    
    output:
    tuple val(segment), path("reference_${segment}.fasta"), emit: reference_segments
    
    script:
    """
    grep -A 1 '${segment}' ${reference} > reference_${segment}.fasta
    """
}