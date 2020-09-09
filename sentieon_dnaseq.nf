#! /usr/bin/env nextflow

fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_alt = file( params.fasta_ref+'.64.alt' )
fasta_ref_amb = file( params.fasta_ref+'.64.amb' )
fasta_ref_ann = file( params.fasta_ref+'.64.ann' )
fasta_ref_bwt = file( params.fasta_ref+'.64.bwt' )
fasta_ref_pac = file( params.fasta_ref+'.64.pac' )
fasta_ref_sa = file( params.fasta_ref+'.64.sa' )
fasta_ref_dict = file( params.fasta_ref.replace(".fasta",".dict") )

fastq1Path = file(params.fastqr1)
fastq2Path = file(params.fastqr2)
sample_id = params.sample_id
threads = params.threads
outpath = params.outpath
code_ver = params.code_ver
iswes = params.wes
aws_region = params.aws_region
jobdef = params.jobdef

reference_mills = file(params.ref_mills)
reference_mills_index = file(params.ref_mills+'.tbi')
reference_dbsnp = file(params.ref_dbsnp)
reference_dbsnp_index = file(params.ref_dbsnp+'.idx')
reference_bedfile = params.wes_bedfile


process alignment {
    tag "${sample_id}"

    cpus 35
    memory '61400 MB'

    input:
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai
    file fasta_ref_alt
    file fasta_ref_amb
    file fasta_ref_ann
    file fasta_ref_bwt
    file fasta_ref_pac
    file fasta_ref_sa
    file fastq1Path
    file fastq2Path

    output:
    file "${sample_id}_sorted.bam" into outputs_sorted_bam
    file "${sample_id}_sorted.bam.bai" into outputs_indexed_bam
    file "code_ver"

    publishDir "${outpath}/fastq/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    bwa mem -v 1 -K '${params.chunk_size}' -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:Illumina" -M -t '${threads}' '${fasta_ref}' '${fastq1Path}' '${fastq2Path}' | sentieon util sort -r '${fasta_ref}' -o '${sample_id}_sorted.bam' -t '${threads}' --sam2bam -i -
    echo ${code_ver} > code_ver
    """
}

process alignment_metrics {
    tag "${sample_id}"

    cpus 18

    input:
    file aligned_bam from outputs_sorted_bam
    file bam_index from outputs_indexed_bam
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai

    output:
    file "${sample_id}.mq_metrics.txt" into output_meanqualitybycycle
    file "${sample_id}.qd_metrics.txt" into output_qualdistribution
    file "${sample_id}.gc_summary.txt" into output_gcsummary
    file "${sample_id}.gc_metrics.txt" into output_gcmetrics
    file "${sample_id}.aln_metrics.txt" into output_alignsummary
    file "${sample_id}.is_metrics.txt" into output_insertsize
    file "code_ver"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon driver --input '${aligned_bam}' \
                                                     --reference '${fasta_ref}'  \
                                                     --thread_count '${threads}' \
                                                     --algo MeanQualityByCycle '${sample_id}.mq_metrics.txt'  \
                                                     --algo QualDistribution '${sample_id}.qd_metrics.txt'  \
                                                     --algo GCBias --summary '${sample_id}.gc_summary.txt' '${sample_id}.gc_metrics.txt'  \
                                                     --algo AlignmentStat '${sample_id}.aln_metrics.txt'  \
                                                     --algo InsertSizeMetricAlgo '${sample_id}.is_metrics.txt'
    echo ${code_ver} > code_ver
    """
}

process plot_alignment_metrics_gc {
    tag "${sample_id}"

    errorStrategy 'finish'

    cpus 18

    input:
    file input_gcmetrics from output_gcmetrics

    output:
    file "${sample_id}.gc_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon plot GCBias -o "${sample_id}.gc_metrics_plot.pdf" '${input_gcmetrics}'
    """
}

process plot_alignment_metrics_qd {
    tag "${sample_id}"

    cpus 18

    input:
    file input_qualdistribution from output_qualdistribution

    output:
    file "${sample_id}.qd_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon plot QualDistribution -o "${sample_id}.qd_metrics_plot.pdf" '${input_qualdistribution}'
    """
}

process plot_alignment_metrics_mq {
    tag "${sample_id}"

    cpus 18

    input:
    file input_meanqualitybycycle from output_meanqualitybycycle

    output:
    file "${sample_id}.mq_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon plot MeanQualityByCycle -o "${sample_id}.mq_metrics_plot.pdf" '${input_meanqualitybycycle}'
    """
}

process plot_alignment_metrics_isize {
    tag "${sample_id}"

    cpus 18

    input:
    file input_insertsize from output_insertsize

    output:
    file "${sample_id}.is_metrics_plot.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon plot InsertSizeMetricAlgo -o "${sample_id}.is_metrics_plot.pdf" '${input_insertsize}'
    """
}

process locus_collector {
    tag "${sample_id}"

    cpus 18

    input:
    file aligned_bam from outputs_sorted_bam
    file bam_index from outputs_indexed_bam

    output:
    file "${sample_id}.score.txt" into output_score_info
    file "${sample_id}.score.txt.idx" into output_score_info_index

    
    """
    sentieon driver --input '${aligned_bam}' \
                                                     --thread_count '${threads}' \
                                                     --algo LocusCollector \
                                                     --fun score_info '${sample_id}.score.txt'

    """
}

process deduplication {
    tag "${sample_id}"

    cpus 18

    input:
    file aligned_bam from outputs_sorted_bam
    file bam_index from outputs_indexed_bam
    file score_info from output_score_info
    file score_info_index from output_score_info_index

    output:
    file "${sample_id}.deduped.bam" into outputs_deduped_bam
    file "${sample_id}.deduped.bam.bai" into outputs_deduped_indexed_bam
    file "${sample_id}.dedup_metrics.txt" into outputs_deduped_metrics

    
    """
    sentieon driver --input '${aligned_bam}' \
                                                     --thread_count '${threads}' \
                                                     --algo Dedup \
                                                     --metrics '${sample_id}.dedup_metrics.txt' \
                                                     --rmdup \
                                                     --score_info '${score_info}' \
                                                     '${sample_id}.deduped.bam'
    """
}

process realignment {
    tag "${sample_id}"

    cpus 18

    input:
    file deduped_bam from outputs_deduped_bam
    file bam_index from outputs_deduped_indexed_bam
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai
    file reference_mills
    file reference_mills_index

    output:
    file "${sample_id}.bam" into outputs_bam_realignment
    file "${sample_id}.bam.bai" into outputs_realigned_indexed_bam
    file "code_ver"

    publishDir "${outpath}/bam/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon driver --input '${deduped_bam}' \
                                                     --reference '${fasta_ref}' \
                                                     --thread_count '${threads}' \
                                                     --algo Realigner \
                                                     --known_sites '${reference_mills}' \
                                                     '${sample_id}.bam'

    echo "${code_ver}" > code_ver
    """
}

process qualcal {
    tag "${sample_id}"

    cpus 11

    input:
    file realigned_bam from outputs_bam_realignment
    file bam_index from outputs_realigned_indexed_bam
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai
    file reference_mills
    file reference_mills_index
    file reference_dbsnp
    file reference_dbsnp_index

    output:
    file "${sample_id}.recal.table" into outputs_recal_table
    
    """
    sentieon driver --input '${realigned_bam}' \
                                                     --reference '${fasta_ref}' \
                                                     --thread_count '${threads}' \
                                                     --algo QualCal \
                                                     --known_sites '${reference_dbsnp}' \
                                                     --known_sites '${reference_mills}' \
                                                     '${sample_id}.recal.table'
    """
}

process qualcalpost {
    tag "${sample_id}"

    cpus 11

    input:
    file realigned_bam from outputs_bam_realignment
    file bam_index from outputs_realigned_indexed_bam
    file qualcal from outputs_recal_table
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai
    file reference_mills
    file reference_mills_index
    file reference_dbsnp
    file reference_dbsnp_index

    output:
    file "${sample_id}.recal.table.post" into outputs_recal_table_post

    """
    sentieon driver --input '${realigned_bam}' \
                                                     --qual_cal '${qualcal}' \
                                                     --reference '${fasta_ref}' \
                                                     --thread_count '${threads}' \
                                                     --algo QualCal \
                                                     --known_sites '${reference_dbsnp}' \
                                                     --known_sites '${reference_mills}' \
                                                     '${sample_id}.recal.table.post'
    """
}

process applyrecal {
    tag "${sample_id}"

    cpus 11

    input:
    file qualcal from outputs_recal_table
    file recal_table from outputs_recal_table_post

    output:
    file "${sample_id}.recal.csv" into outputs_bqsr
    
    """
    sentieon driver --thread_count '${threads}' \
                                                     --algo QualCal \
                                                     --after '${recal_table}' \
                                                     --before '${qualcal}' \
                                                     --plot '${sample_id}.recal.csv'
    """
}

process plotbqsr {
    tag "${sample_id}"
    
    cpus 10

    input:
    file recal_table from outputs_bqsr

    output:
    file "${sample_id}.recal.pdf"

    publishDir "${outpath}/metrics/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    sentieon plot QualCal -o '${sample_id}.recal.pdf' '${recal_table}'
    """
}

process haplotypecaller {
    tag "${sample_id}"

    cpus 30

    input:
    file realigned_bam from outputs_bam_realignment
    file bam_index from outputs_realigned_indexed_bam
    file recal_table from outputs_recal_table
    file fasta_ref
    file fasta_ref_dict
    file fasta_ref_fai
    file reference_dbsnp
    file reference_dbsnp_index

    output:
    file "${sample_id}.g.vcf.gz" into outputs_gvcf
    file "${sample_id}.g.vcf.gz.tbi" into outputs_gvcf_index
    file "code_ver"

    publishDir "${outpath}/gvcf/${sample_id}/", mode: 'copy', overwrite: false
    
    """
    ARG_INTERVAL=""
    if [[ ${iswes} == "true" ]]
    then
        aws s3 cp ${reference_bedfile} wes.bed --region ${aws_region}
        ARG_INTERVAL="--interval wes.bed"
    fi
    sentieon driver \$ARG_INTERVAL --input '${realigned_bam}' \
                                                               --qual_cal '${recal_table}' \
                                                               --reference '${fasta_ref}' \
                                                               --thread_count '${threads}' \
                                                               --algo Haplotyper \
                                                               --call_conf 30 \
                                                               --dbsnp '${reference_dbsnp}' \
                                                               --emit_conf 30 \
                                                               --emit_mode Gvcf \
                                                               '${sample_id}.g.vcf.gz'

    echo "${code_ver}" > code_ver
    """
    
}
