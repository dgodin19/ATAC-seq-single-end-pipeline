include {GET_ACCESSIONS} from './modules/get_accessions'
include {DOWNLOAD_SAMPLES} from './modules/download_data'
include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {BOWTIE2_BUILD} from './modules/build'
include {BOWTIE2_ALIGN} from './modules/align'
include {BAMCOVERAGE} from './modules/deeptools_bamcoverage'
include {COMPUTEMATRIX as COMPUTEMATRIX_CDC1} from './modules/deeptools_computematrix'
include {COMPUTEMATRIX as COMPUTEMATRIX_CDC2} from './modules/deeptools_computematrix'
include {PLOTPROFILE as PLOTPROFILE_CDC1} from './modules/deeptools_plotprofile'
include {PLOTPROFILE as PLOTPROFILE_CDC2} from './modules/deeptools_plotprofile'
include {COMPUTE_TSS_ENRICH as COMPUTE_TSS_ENRICH_CDC1} from './modules/compute_tss_enrich'
include {COMPUTE_TSS_ENRICH as COMPUTE_TSS_ENRICH_CDC2} from './modules/compute_tss_enrich'
include {PLOTHEATMAP as PLOTHEATMAP_CDC1} from './modules/deeptools_plotheatmap'
include {PLOTHEATMAP as PLOTHEATMAP_CDC2} from './modules/deeptools_plotheatmap'
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

    // Map to include celltype
    diffbind_gain_loss_ch = DIFFBIND.out.gain_loss
    .map { gain, loss ->
        def m = gain.baseName =~ /(cDC\d+)/
        def celltype = m ? m[0][1] : 'UNKNOWN'
        tuple(celltype, gain, loss)
    }
    

    run2cond = [
    'SRR28895183': 'cDC2_KO_2',
    'SRR28895184': 'cDC1_KO_2',
    'SRR28895185': 'cDC2_WT_2',
    'SRR28895186': 'cDC1_WT_2',
    'SRR28895187': 'cDC2_KO_1',
    'SRR28895188': 'cDC1_KO_1',
    'SRR28895189': 'cDC2_WT_1',
    'SRR28895190': 'cDC1_WT_1'
    ]

    BAMCOVERAGE(BOWTIE2_ALIGN.out.filtered_bam, run2cond)

    branch_ch = BAMCOVERAGE.out.branch { it ->
    cDC1_WT: it[1].name.contains('cDC1_WT')
    cDC1_KO: it[1].name.contains('cDC1_KO')
    cDC2_WT: it[1].name.contains('cDC2_WT')
    cDC2_KO: it[1].name.contains('cDC2_KO')
    }
   
    sigpeak_ch = DIFFBIND.out.significant_peaks.map { path ->
    def m = path.baseName =~ /(cDC\d+)/
    def celltype = m ? m[0][1] : 'UNKNOWN'
    tuple(celltype, path)
    }
    sigpeak_ch.view()


    // Collect bigwig lists
    branch_ch.cDC1_WT.collect() | set { cDC1_WT_bw }
    branch_ch.cDC1_KO.collect() | set { cDC1_KO_bw }
    branch_ch.cDC2_WT.collect() | set { cDC2_WT_bw }
    branch_ch.cDC2_KO.collect() | set { cDC2_KO_bw }

    // Merge WT + KO per cell type
    cDC1_all_bw = cDC1_WT_bw.combine(cDC1_KO_bw).map { it.flatten() }
    cDC2_all_bw = cDC2_WT_bw.combine(cDC2_KO_bw).map { it.flatten() }
    

    cDC1_all_bw_redux = cDC1_all_bw.map { list -> list.collate(2).collect { it[1] } }
    cDC2_all_bw_redux = cDC2_all_bw.map { list -> list.collate(2).collect { it[1] } }
    
    cDC1_sigpeak = sigpeak_ch.filter { it[0] == 'cDC1' }.map { it[1] }
    cDC2_sigpeak = sigpeak_ch.filter { it[0] == 'cDC2' }.map { it[1] }

    COMPUTEMATRIX_CDC1('cDC1',cDC1_all_bw_redux, cDC1_sigpeak)
    COMPUTEMATRIX_CDC2('cDC2', cDC2_all_bw_redux, cDC2_sigpeak)

    PLOTPROFILE_CDC1(COMPUTEMATRIX_CDC1.out)
    PLOTPROFILE_CDC2(COMPUTEMATRIX_CDC2.out)

    PLOTHEATMAP_CDC1(COMPUTEMATRIX_CDC1.out)
    PLOTHEATMAP_CDC2(COMPUTEMATRIX_CDC2.out)

    COMPUTE_TSS_ENRICH_CDC1(PLOTPROFILE_CDC1.out.signal_tab)
    COMPUTE_TSS_ENRICH_CDC2(PLOTPROFILE_CDC2.out.signal_tab)

    
    genome_idx_ch = SAMTOOLS_GENOME(params.genome)
    ANNOTATE(DIFFBIND.out.significant_peaks, genome_idx_ch, params.gtf)
    FIND_MOTIFS(DIFFBIND.out.significant_peaks, genome_idx_ch)

    
}