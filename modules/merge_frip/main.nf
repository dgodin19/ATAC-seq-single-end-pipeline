
#!/usr/bin/env nextflow

process MERGE_FRIP {
    label 'process_medium'
    publishDir params.outdir, mode: "copy", pattern: '*.txt'
    
    input:
    path frip_files

    output:
    path("frip_summary.txt")

    script:
    """
    echo "run,frip" > frip_summary.txt
    for f in ${frip_files}; do
        tail -n +2 \$f >> frip_summary.txt
    done
    """
}
