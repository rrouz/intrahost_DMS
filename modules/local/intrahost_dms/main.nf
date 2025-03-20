process INTRAHOST_DMS {
    publishDir "${params.outdir}/intrahost_dms", mode: 'copy'

    input:
    path dms_files

    output:
    path "combined_variants.tsv", emit: combined_variants

    script:
    """
    intrahost_dms.py \
        --dms_file 'https://raw.githubusercontent.com/dms-vep/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/main/results/summaries/phenotypes.csv' \
        --data_dir ${projectDir}/data/variants \
        --output_dir ./
    """
}