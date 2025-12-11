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
include {PREPARE_DIFFBIND_SAMPLESHEET} from './modules/prepare_diffbind_samplesheet'
include {MERGE_SAMPLESHEETS} from './modules/merge_samplesheets'
include {DIFFBIND} from './modules/diffbind'
/*include {MERGE_PEAKS as MERGE_PEAKS_KO} from './modules/merge_peaks'
include {MERGE_PEAKS as MERGE_PEAKS_WT} from './modules/merge_peaks'
include {ANNOTATE as ANNOTATE_KO} from './modules/annotate'
include {ANNOTATE as ANNOTATE_WT} from './modules/annotate'
include {FIND_MOTIFS as FIND_MOTIFS_KO} from './modules/find_motifs'
include {FIND_MOTIFS as FIND_MOTIFS_WT} from './modules/find_motifs'*/


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
    prepared_ch.view()
    
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

    /*merged_csv_ch.branch {
    cdc1: it.name.contains('cDC1')
    cdc2: it.name.contains('cDC2')
    }.set { branched }

    DIFFBIND_CDC1(branched.cdc1)
    DIFFBIND_CDC2(branched.cdc2)*/

    // Output from callpeaks
    //tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),
    //      path("${run}_peaks.narrowPeak"),
    //      path("${run}_summits.bed"),
    //      path("${run}_peaks.xls")

    //Output from align
    //tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), 
    //path("${run}.filtered.bam"), path("${run}.filtered.bam.bai"), emit: filtered_bam

    //Structure of metadata
    // run,genotype,library_name,condition

    //Target structure of input csv file for cDC1
    // SampleID	Tissue	Factor	Condition	Replicate	bamReads	Peaks	PeakCaller
    // SRR28895184	cells	ATAC	cDC1_KO	2	SRR28895184.filtered.bam	SRR28895184_peaks.narrowPeak	narrowPeak
    // SRR28895186	cells	ATAC	cDC1_WT	2	SRR28895186.filtered.bam	SRR28895186_peaks.narrowPeak	narrowPeak
    // SRR28895188	cells	ATAC	cDC1_KO	1	SRR28895188.filtered.bam	SRR28895188_peaks.narrowPeak	narrowPeak
    //SRR28895190	cells	ATAC	cDC1_WT	1	SRR28895190.filtered.bam	SRR28895190_peaks.narrowPeak	narrowPeak

    //Target structure of input csv file for cDC2
    //SampleID	Tissue	Factor	Condition	Replicate	bamReads	Peaks	PeakCaller
    // SRR28895183	cells	ATAC	cDC2_KO	2	SRR28895183.filtered.bam	SRR28895183_peaks.narrowPeak	narrowPeak
    // SRR28895185	cells	ATAC	cDC2_WT	2	SRR28895185.filtered.bam	SRR28895185_peaks.narrowPeak	narrowPeak
    // SRR28895187	cells	ATAC	cDC2_KO	1	SRR28895187.filtered.bam	SRR28895187_peaks.narrowPeak	narrowPeak
    // SRR28895189	cells	ATAC	cDC2_WT	1	SRR28895189.filtered.bam	SRR28895189_peaks.narrowPeak	narrowPeak

    
}