#!/usr/bin/env nextflow

process GET_ACCESSIONS {
    label 'process_single'
    conda 'envs/ncbidatasets_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.csv'

    input:
    val project_accession

    output:
    path 'accessions.csv'

    script:
    """
    esearch -db sra -query ${project_accession} | efetch -format runinfo > accessions.csv
    """

    stub:
    """
    echo "Run,LibraryLayout,LibrarySource,cell_type,genotype" > accessions.csv
    echo "SRR123456,SINGLE,GENOMIC,cell_type_example,wildtype" >> accessions.csv
    """
}