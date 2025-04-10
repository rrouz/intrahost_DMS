process GOFASTA_VARIANTS {
    publishDir "${params.outdir}/mutations", mode: 'copy'
    
    input:
    tuple val(segment), path(reference_fasta), path(combined_fasta)
    
    output:
    tuple val(segment), path("aa_changes_${segment}.csv"), emit: aa_changes
    
    script:
    """
    minimap2 -t ${task.cpus} \\
             -a \\
             -x asm20 \\
             --score-N=0 \\
             -N 0 \\
             -p 0.8 \\
             ${reference_fasta} \\
             ${combined_fasta} > aligned_${segment}.sam
             
    if [ \$? -ne 0 ]; then
        echo "Error: minimap2 alignment failed"
        exit 1
    fi
    
    # Ensure gofasta is in PATH and correctly running with your fork's options
    which gofasta
    gofasta -h || true
    
    gofasta sam variants \\
        -s aligned_${segment}.sam \\
        --annotation "${params.data_dir}/gff/${segment}.gff" \\
        --append-snps \\
        --append-codons \\
        -o aa_changes_${segment}.csv
        
    if [ \$? -ne 0 ]; then
        echo "Error: gofasta variants analysis failed"
        exit 1
    fi
    """
}