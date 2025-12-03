#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_low'
	conda 'envs/multiqc_env.yml'
	publishDir params.outdir, mode: "copy", pattern: '*.html'
	
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