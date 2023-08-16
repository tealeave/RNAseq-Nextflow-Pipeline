#!/usr/bin/env nextflow

WORK_DIR="/dfs6/pub/ddlin/kaiser_lab/my_nextflow/downsampled_fastq"
GENOME="/dfs6/pub/ddlin/refs/hg38.fa"
GTF="/dfs6/pub/ddlin/refs/hg38_gencode.v44.annotation.gtf"
// STAR_INDEX="/dfs6/pub/ddlin/refs/hg38_star_index/"
// Generate both STAR and RSEM index
// rsem-prepare-reference -p 32 --gtf hg38_gencode.v44.annotation.gtf --star hg38.fa hg38_RSEM_ref/hg38_RSEM_ref 
STAR_INDEX="/dfs6/pub/ddlin/refs/hg38_RSEM_ref/"
RSEM_REF="/dfs6/pub/ddlin/refs/hg38_RSEM_ref/hg38_RSEM_ref"
RMATS="/dfs6/pub/ddlin/kaiser_lab/my_nextflow/rmats-turbo/rmats.py"
RMATS_ENV="/pub/ddlin/kaiser_lab/my_nextflow/rmats-turbo/conda_envs/rmats"
RMATS_BAM_TXT1="/dfs6/pub/ddlin/kaiser_lab/my_nextflow/MB468_000_REP1_bam.txt"
RMATS_BAM_TXT2="/dfs6/pub/ddlin/kaiser_lab/my_nextflow/R8_000_REP1_bam.txt"

nextflow.enable.dsl=2

process trimGalore {
    publishDir "trimmed_reads", mode: 'copy'

    input:
    tuple val(sample), path(read)

    output:
    tuple val(sample), path("${sample}_trimmed.fq.gz"), emit: trimmed

    script:
    """
    trim_galore --fastqc --output_dir . $read
    """
}


process star {
    publishDir "star_aligned", mode: 'copy'

    input:
    tuple val(sample), path(read)

    output:
    tuple val(sample), path("${sample}_Aligned.toTranscriptome.out.bam"), path("${sample}_Aligned.sortedByCoord.out.bam") , emit: aligned

    script:
    """
    STAR --readFilesCommand zcat \
         --genomeDir $STAR_INDEX \
         --readFilesIn $read \
         --runThreadN 12 \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode TranscriptomeSAM \
         --outFileNamePrefix ${sample}_
    """
}


process rsem {
    publishDir "rsem_output", mode: 'copy'

    input:
    tuple val(sample), path(transcriptome_alignment), path(genome_alignment)

    output:
    tuple path("${sample}_quantified.genes.results"), path("${sample}_quantified.isoforms.results"), emit: quantified

    script:
    """
    rsem-calculate-expression \\
        --num-threads 12 \\
        --forward-prob 0.5 \\
        --bam $transcriptome_alignment \\
        --no-bam-output \\
        $RSEM_REF \\
        ${sample}_quantified
    """
}

process rmats {
    publishDir "rmats_output", mode: 'copy'

    input:
    tuple path(bam_txt1), path(bam_txt2)

    output:
    path("${bam_txt1.baseName}_vs_${bam_txt2.baseName}_rmats_output"), emit: rmats_out

    script:
    """
    conda activate $RMATS_ENV
    python  $RMATS\\
        --b1 $bam_txt1 \\
        --b2 $bam_txt2 \\
        -t single \\
        --readLength 100 \\
        --gtf $GTF \\
        --od ${bam_txt1.baseName}_vs_${bam_txt2.baseName}_rmats_output \\
        --nthread 18
    """
}


Channel.fromFilePairs("$WORK_DIR/*.fq.gz", size: 1)
       .map { sample, read -> [sample.replaceFirst(/\.fq\.gz$/, ''), read] }
       .set { reads_ch }

// Define a channel for BAM text file pairs
Channel.fromPath([RMATS_BAM_TXT1, RMATS_BAM_TXT2])
       .set { bam_txt_ch }

workflow {
    reads_ch
        .set { input_ch }
    
    trimGalore(input_ch)
        .set { trimmed_ch }

    star(trimmed_ch)
        .set { star_ch }

    rsem(star_ch)

    // rmats(bam_txt_ch)
}
