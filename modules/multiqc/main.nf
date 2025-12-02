#!/usr/bin/env nextflow

process MULTIQC {
    container 'ghcr.io/bf528/multiqc:latest'
	publishDir params.outdir, mode: "copy", pattern: '*.html'
	label 'process_low'

	input:
	path ('*')

	output:
	path('multiqc_report.html')

	script:
	"""
	multiqc . 
	"""
    stub:
    """
    touch multiqc_report.html
    """
}