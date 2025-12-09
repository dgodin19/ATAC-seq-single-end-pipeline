include {GET_ACCESSIONS} from './modules/get_accessions'
include {DOWNLOAD_SAMPLES} from './modules/download_data'
include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {BOWTIE2_BUILD} from './modules/build'
include {BOWTIE2_ALIGN} from './modules/align'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {MULTIQC} from './modules/multiqc'
include {REMOVE_BLACKLIST} from './modules/remove_blacklist'
include {CALL_PEAKS} from './modules/call_peaks'
include {MERGE_PEAKS as MERGE_PEAKS_KO} from './modules/merge_peaks'
include {MERGE_PEAKS as MERGE_PEAKS_WT} from './modules/merge_peaks'


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
    REMOVE_BLACKLIST_OUT = REMOVE_BLACKLIST(CALL_PEAKS.out, params.blacklist)

    // Create metadata channel
    metadata = Channel.fromPath(params.metadata)
    .splitCsv(header: true)
    .map { row -> 
        [
            row.run, 
            row.condition.contains('KO') ? 'KO' : 'WT'
        ] 
    }

    // Extract bed files from REMOVE_BLACKLIST output
    bed_files = REMOVE_BLACKLIST_OUT.map { tuple -> 
    [tuple[0], tuple[6]]  // Include run ID with bed file
    }

    // Join bed files with metadata
    bed_files_with_condition = bed_files.join(metadata)
    .map { run, bed_file, condition -> 
        [bed_file, condition] 
    }

    // Split bed files by condition
    knockout_bed_files = bed_files_with_condition
    .filter { it[1] == 'KO' }
    .map { it[0] }

    wildtype_bed_files = bed_files_with_condition
    .filter { it[1] == 'WT' }
    .map { it[0] }

    // Collect bed files for each condition
    knockout_bed_files_grouped = knockout_bed_files.collect()
    wildtype_bed_files_grouped = wildtype_bed_files.collect()

    // Merge peaks for each condition
    MERGE_PEAKS_KO(knockout_bed_files_grouped, Channel.value('KO'))
    MERGE_PEAKS_WT(wildtype_bed_files_grouped, Channel.value('WT'))

}