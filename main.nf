include {GET_ACCESSIONS} from './modules/get_accessions'
include {DOWNLOAD_SAMPLES} from './modules/download_data'
include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {BOWTIE2_BUILD} from './modules/build'
include {BOWTIE2_ALIGN} from './modules/align'
include {BAMCOVERAGE} from './modules/deeptools_bamcoverage'
include {COMPUTEMATRIX} from './modules/deeptools_computematrix'
include {PLOTPROFILE} from './modules/deeptools_plotprofile'
include {COMPUTE_TSS_ENRICH} from './modules/compute_tss_enrich'
include {FRIP_CALC} from './modules/frip'
include {MERGE_FRIP} from './modules/merge_frip'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {MULTIQC} from './modules/multiqc'
include {REMOVE_BLACKLIST} from './modules/remove_blacklist'
include {CALL_PEAKS} from './modules/call_peaks'
include {PREPARE_DIFFBIND_SAMPLESHEET} from './modules/prepare_diffbind_samplesheet'
include {MERGE_SAMPLESHEETS} from './modules/merge_samplesheets'
include {DIFFBIND} from './modules/diffbind'
include {SAMTOOLS_GENOME} from './modules/samtools_genome'
include {ANNOTATE} from './modules/annotate'
include {FIND_MOTIFS} from './modules/find_motifs'


workflow {
    // Get accessions with full metadata
    accessions_ch = GET_ACCESSIONS(params.project_accession)
    
    // Convert stdout to file
    accessions_ch
    .splitCsv(header:true)
    .map{ row -> tuple(row.Run, row.BioSample, row.SampleName, row.LibraryLayout, row.LibrarySource, row.Experiment) }
    .set { read_ch }
    DOWNLOAD_SAMPLES(read_ch)
    FASTQC(DOWNLOAD_SAMPLES.out)
    TRIM(DOWNLOAD_SAMPLES.out, params.adapters)
    BOWTIE2_BUILD(tuple("mm10", params.genome))
    BOWTIE2_ALIGN(TRIM.out.trim, BOWTIE2_BUILD.out.index)

    SAMTOOLS_FLAGSTAT(BOWTIE2_ALIGN.out.samtools_flagstat)
    

    multiqc_ch = Channel
    .empty()
    .mix(TRIM.out.log.map { it[6] })  
    .mix(SAMTOOLS_FLAGSTAT.out.flagstat.map { it[6] }) 
    .mix(FASTQC.out.html.map { it[6] })
    .mix(FASTQC.out.zip.map { it[6] })  // Add FastQC zip files
    .flatten()
    .unique()
    .map { 
        println "MultiQC Input: $it"  // Debug print
        it 
    }
    .collect()

    MULTIQC(multiqc_ch)

    CALL_PEAKS(BOWTIE2_ALIGN.out.filtered_bam)
    BAMCOVERAGE(BOWTIE2_ALIGN.out.filtered_bam)
    bigwig_summary = BAMCOVERAGE.out.map { it[1] }.collect()
    COMPUTEMATRIX(bigwig_summary, params.ucsc_genes, params.window)
    PLOTPROFILE(COMPUTEMATRIX.out)
    COMPUTE_TSS_ENRICH(PLOTPROFILE.out.signal_tab)
    
    // Load metadata
    metadata_ch = Channel
        .fromPath(params.metadata)
        .splitCsv(header: true)
        .map { row ->
            def parts = row.condition.split('_')
            tuple(row.run, parts[0], parts[1], parts[2], parts[3])
        }
    

    align_peaks_ch = BOWTIE2_ALIGN.out.filtered_bam
        .join(CALL_PEAKS.out, by: [0,1,2,3,4,5])
    
    FRIP_CALC(align_peaks_ch)
    MERGE_FRIP( FRIP_CALC.out.collect() )

    
    combined_ch = align_peaks_ch.join(metadata_ch, by: 0)
    
    

    prepared_ch = combined_ch.map {
    tuple(
        it[0],   // run
        it[1],   // biosample
        it[2],   // samplename
        it[3],   // library_layout
        it[4],   // library_source
        it[5],   // experiment
        it[6],   // filtered_bam (.bam)
        it[7],   // filtered_bam_index (.bam.bai)
        it[8],   // narrowPeak (.narrowPeak)
        it[11],  // factor ("ATAC")
        it[12],  // cell_type (e.g., "cDC2")
        it[13],  // ko_wt (e.g., "KO")
        it[14]   // replicate (e.g., "1")
    )
    }
    
    
    sheet_ch = PREPARE_DIFFBIND_SAMPLESHEET(prepared_ch)
    
    sheet_ch
    .map { file -> 
    def match = file.name =~ /(cDC\d)/
    def cell_type = match ? match[0][1] : 'UNKNOWN'
    tuple(cell_type, file)
    }
    .groupTuple()
    .set { merged_groups }

    merged_csv_ch = MERGE_SAMPLESHEETS(merged_groups)

    merged_csv_ch
    .map { it }   // identity
    .set { all_sheets }

    DIFFBIND(all_sheets)
    genome_idx_ch = SAMTOOLS_GENOME(params.genome)
    ANNOTATE(DIFFBIND.out, genome_idx_ch, params.gtf)
    FIND_MOTIFS(DIFFBIND.out, genome_idx_ch)

    
}