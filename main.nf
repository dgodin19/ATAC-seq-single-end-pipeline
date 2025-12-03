include {GET_ACCESSIONS} from './modules/get_accessions'
include {DOWNLOAD_SAMPLES} from './modules/download_data'
include {FASTQC} from './modules/fastqc'
include {TRIM} from './modules/trimmomatic'
include {BOWTIE2_BUILD} from './modules/build'
include {BOWTIE2_ALIGN} from './modules/align'
include {SAMTOOLS_FLAGSTAT} from './modules/samtools_flagstat'
include {SHIFT_ALIGN} from './modules/shift_align'

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
    SHIFT_ALIGN(BOWTIE2_ALIGN.out.filtered_bam)


    


}