process PULL_DMS {
    publishDir "${params.outdir}/dms", mode: 'copy'
    
    input:
    tuple val(segment), path(aa_changes_csv)
    
    output:
    tuple val(segment), path("dms_${segment}.json"), path("dms_${segment}.csv"), emit: dms_data
    
    script:
    """
    pull_dms.py \
        --dms-file 'https://raw.githubusercontent.com/dms-vep/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/main/results/summaries/phenotypes.csv' \
        --dms-file-2 'https://raw.githubusercontent.com/dms-vep/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/main/results/summaries/phenotypes_per_antibody_escape.csv' \
        --h5_site_header "sequential_site" \
        --h3_site_header "site" \
        --mutation_file ${aa_changes_csv} \
        --output_file dms_${segment}.json
    """
}