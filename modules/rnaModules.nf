#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

inhouse_genelist="/data/shared/genomes/databases/genelists/tumortarget/240123.inhouse.MOMA.Fusion.241genes.for.grep.txt"
inhouse_splicing_genelist="/data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes_for_grep.txt"

//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    default:
    gatk_image="gatk4400.sif";
    break;
}

switch (params.server) {
    case 'lnx01':
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/data/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        refFilesDir="/data/shared/genomes";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
    break;
    default:
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/,/lnx01_data4/:/lnx01_data4/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        refFilesDir="/fast/shared/genomes";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml"
    break;
}

switch (params.genome) {
    case 'hg19':
        // Genome assembly files:
        genome_fasta = "${refFilesDir}/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "${refFilesDir}/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "${refFilesDir}/hg19/human_g1k_v37.dict"

        // Gene and transcript annotation files:
        gene_bed12 = "${refFilesDir}/hg19/gene.annotations/hg19.refGene.201228.BED12.bed"
        transcript_fasta="${refFilesDir}/GRCh37/gene.annotations/gencode.v19.pc_transcripts.fa"
        //gencode_gtf = "${refFilesDir}/GRCh37/gene.annotations/gencode.v19.annotation.gtf"
        gencode_gtf = "${refFilesDir}/hg19/gene.annotations/hg19.refGene.201228.gtf" // make sure to change this after testing!
         //Program  files:
        msisensor_list="${refFilesDir}/hg19/human_g1k_v37.microsatellites.list"
        genome_lib_starfusion= "${refFilesDir}/hg19/CTAT/GRCh37.Apr032020.PNP/ctat_genome_lib_build_dir/"
        arriba_blacklist = "/data/shared/programmer/arriba_v2.0.0/database/blacklist_hg19_hs37d5_GRCh37_v2.0.0.tsv.gz"
        arriba_cytoband= "/data/shared/genomes/hg19/arriba/cytobands_hg19_hs37d5_GRCh37_2018-02-23.tsv"
        arriba_protein_gff = "/data/shared/genomes/hg19/arriba/protein_domains_hg19_hs37d5_GRCh37_2019-07-05.gff3"

        // Program indexes
        index_rsem = "/data/shared/genomes/GRCh37/rsem/rsem_gencode"
        //index_star = "/data/shared/genomes/GRCh37/star_align_grch37"
        index_star = "/data/shared/genomes/hg19/star_align_hg19"

        //regions:
        qualimap_ROI="/data/shared/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.merged.sorted.6col.bed"
        ROI="/data/shared/genomes/hg19/interval.files/WES/IDT.exomes.EV7/EV7.ROI.bed"

        break;
    case 'hg38':
        // Genome assembly files:
        genome_fasta = "${refFilesDir}/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38.primary.dict"

        // Gene and transcript annotation files:
        gene_bed12="${refFilesDir}/hg38/gene.annotations/gencode.v36.BED12.bed"
        transcript_fasta="${refFilesDir}/hg38/gene.annotations/gencode.v36.transcripts.fa"
        gencode_gtf = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.gff3"
        gencode_gtf_collapsed="${refFilesDir}/hg38/gene.annotations/gencode.v36.annotation.collapsed.gtf"
        ensemble102_gtf="${refFilesDir}/hg38/gene.annotations/ensembl102/Homo_sapiens.GRCh38.102.gtf"
        ensemble102_transcript_fasta="${refFilesDir}/hg38/gene.annotations/ensembl102/Homo_sapiens.GRCh38.102.cdna.all.fa.gz"
        //Program  files:
        hmftools_data_dir="/data/shared/genomes/hg38/program_DBs/hmftools/hmf_pipeline_resources.38_v6.0/"


        msisensor_list="/data/shared/genomes/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        genome_lib_starfusion="/data/shared/genomes/hg38/program_DBs/CTAT/GRCh38_v33_Apr062020.PNP/ctat_genome_lib_build_dir/"
        arriba_cytoband="/data/shared/genomes/hg38/program_DBs/arriba/cytobands_hg38_GRCh38_v2.1.0.tsv"
        arriba_blacklist = "/data/shared/genomes/hg38/program_DBs/arriba/blacklist_hg38_GRCh38_v2.1.0.tsv.gz"
        arriba_protein_gff="/data/shared/genomes/hg38/program_DBs/arriba/protein_domains_hg38_GRCh38_v2.1.0.gff3"
        arriba_known_fusions="/data/shared/genomes/hg38/program_DBs/arriba/known_fusions_hg38_GRCh38_v2.1.0.tsv.gz"

        //fusioncallers
        fusioncatcher_db="/data/shared/genomes/hg38/program_DBs/fusioncatcher/human_v102/"
        fusionreport_db="/data/shared/genomes/hg38/program_DBs/fusion_report_db/"
        // Program indexes:
        index_rsem = "/data/shared/genomes/hg38/rsem/rsem_hg38_gencode36"
        index_star = "${refFilesDir}/hg38/STAR/"
        kallisto_index="/data/shared/genomes/hg38/program_DBs/kallisto/Homo_sapiens.GRCh38.102.cdna.all_kallisto_K31.idx"
        //regions:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.6col.bed"
        ROI="/data/shared/genomes/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        gencode36_coding_exons="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"
        break;
}

outputDir="Results"

/*
Default input channel
meta: caseID, npnNormal,npnTumorDNA,npnTumorRNA,pcgr
data (default): 0:cramNormal, 1:craiNormal, 2:cramTumorDNA, 3:craiTumorDNA, 4:R1_RNA, 5:R2_RNA

NEW DATA STRUCTURE:
data rnaInput (r1, r2)

*/


process inputLinks {
    errorStrategy 'ignore'
    publishDir "${meta.id}/inputSymlinksFASTQ/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path(data)
    script:
    """
    """
}

process fastp_TRIM {
    publishDir "${meta.id}/toolsOutput/QC/", mode: 'copy', pattern: '*.{html,json}'
    cpus 10
    tag "$meta.id"

    input:
    tuple val(meta), path(data)
    output:
    path("*.{html,json}")                                           ,emit: fastp_results
    tuple val(meta), path("*.fastp.fq.gz")    ,emit: fastp_out_ch

    script:
    """
    singularity run -B ${s_bind} \
    ${simgpath}/fastp.sif \
    -i ${data[0]} -I ${data[1]} \
    -o ${data[0].baseName}.fastp.fq.gz -O ${data[1].baseName}.fastp.fq.gz \
    --json ${meta.npnTumorRNA}.fastp.json \
    --html ${meta.npnTumorRNA}.fastp.html \
    -w ${task.cpus}
    """
}

process align_STAR {
    tag "$meta.id"
    cpus 40
    maxForks 3
    publishDir "${meta.id}/alignments/", mode: 'copy', pattern: '*.hg38.cra*'
    publishDir "${meta.id}/toolsOutput/QC/", mode: 'copy', pattern: '*.Log.*'
 
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/star2711b'

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.STAR.*.cram"),path("${meta.id}.${meta.npnTumorRNA}_RNA.STAR.*.crai"), emit: star_out_cram_ch 
    tuple val(meta), path("*.hg38.cra*")
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.STAR.Aligned.toTranscriptome.*.bam"), emit:rsem_input_bam    // into rsem_input_bam

    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.STAR_chimeric.hg38.cram"),path("${meta.id}.${meta.npnTumorRNA}_RNA.STAR_chimeric.*.crai"), emit: isofox_input_cram 

    tuple val(meta), path("*.Chimeric.out.junction"),emit: chimeric_junctions_out  

    //tuple val(meta), path("*.STAR.SJ.out.tab"),emit: star_sjtab_out
    tuple val(meta), path("*.STAR.SJ.out.tab"), path("*.Chimeric.out.junction"), emit: sj_chimJunction    
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.forArriba.*.bam"),path("${meta.id}.${meta.npnTumorRNA}_RNA.forArriba.*.bai"), emit: arriba_input_bam

    path("*.STAR.Log.*")
    path("*.forArriba.Log.out"), emit: trinity_collect // into (trinity_collect_ch1,trinity_collect_ch2)
    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${index_star} \
        --sjdbGTFfile ${gencode_gtf} \
        --readFilesIn ${data} \
        --outFileNamePrefix ${meta.id}.${meta.npnTumorRNA}_RNA.STAR. \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif --outSAMunmapped Within \
        --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${meta.id}.${meta.npnTumorRNA}_RNA \
        --chimSegmentMin 12 \
        --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
        --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
        --outReadsUnmapped None \
        --alignInsertionFlush Right \
        --alignSplicedMateMapLminOverLmate 0 \
        --alignSplicedMateMapLmin 30 \
        --quantMode TranscriptomeSAM GeneCounts \
        --readFilesCommand zcat
 
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${index_star} \
        --sjdbGTFfile ${gencode_gtf} \
        --readFilesIn ${data} \
        --outFileNamePrefix ${meta.id}.${meta.npnTumorRNA}_RNA.forArriba. \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMstrandField intronMotif --outSAMunmapped Within \
        --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${meta.id}.${meta.npnTumorRNA}_RNA \
        --chimSegmentMin 12 \
        --chimOutType WithinBAM \
        --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG -4 \
        --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 \
        --chimJunctionOverhangMin 12 --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
        --outReadsUnmapped None \
        --readFilesCommand zcat
    
    bamtools index -in ${meta.id}.${meta.npnTumorRNA}_RNA.forArriba.Aligned.sortedByCoord.out.bam
    bamtools index -in ${meta.id}.${meta.npnTumorRNA}_RNA.STAR.Aligned.sortedByCoord.out.bam
    
    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${meta.npnTumorRNA}_RNA.STAR.hg38.cram \
    ${meta.id}.${meta.npnTumorRNA}_RNA.STAR.Aligned.sortedByCoord.out.bam

    samtools index ${meta.id}.${meta.npnTumorRNA}_RNA.STAR.hg38.cram

    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${meta.npnTumorRNA}_RNA.STAR_chimeric.hg38.cram \
    ${meta.id}.${meta.npnTumorRNA}_RNA.forArriba.Aligned.sortedByCoord.out.bam

    samtools index ${meta.id}.${meta.npnTumorRNA}_RNA.STAR_chimeric.hg38.cram 
    """
}

process align_STAR_hmf {
    tag "$meta.id"
    cpus 40
    maxForks 3
//    publishDir "${meta.id}/CRAM/Isofox/", mode: 'copy', pattern: '*.hg38.cra*'
//    publishDir "${meta.id}/CRAM/Isofox/", mode: 'copy', pattern: '*.hg38.ba*'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/star2711b'

    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.hg38.cram"), path("${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.*.crai"), emit: isofox_input_cram 
    path("*.out.ba*")

    script:
    """
    STAR --runThreadN ${task.cpus} \
        --genomeDir ${index_star} \
        --sjdbGTFfile ${gencode_gtf} \
        --readFilesIn ${data} \
        --outFileNamePrefix ${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox. \
        --outSAMattrRGline ID:${date}.${user} LB:library PL:illumina PU:illumina SM:${meta.id}.${meta.npnTumorRNA}_RNA \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes All \
        --outFilterMultimapNmax 10 \
        --outFilterMismatchNmax 3 \
        --limitOutSJcollapsed 3000000 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM SoftClip \
        --chimJunctionOverhangMin 10 \
        --chimSegmentReadGapMax 3 \
        --chimScoreMin 1 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --outFilterMatchNmin 35 \
        --alignSplicedMateMapLminOverLmate 0.33 \
        --alignSplicedMateMapLmin 35 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outReadsUnmapped None \
        --readFilesCommand zcat
    
    bamtools index -in ${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.Aligned.sortedByCoord.out.bam

    samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.hg38.cram \
    ${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.Aligned.sortedByCoord.out.bam

    samtools index ${meta.id}.${meta.npnTumorRNA}_RNA.forIsofox.hg38.cram
    """
}

process isofox {
    tag "$meta.id"
    cpus 10
    maxForks 3
    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    publishDir "${meta.id}/toolsOutput/genecounts/", mode: 'copy'
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/isofox171'

    input:
    tuple val(meta), path(data) //data: [cram,crai]
    output:
    tuple val(meta),path("isofox/"), emit: isofoxDir
    script:
    """
    mkdir isofox/
    isofox "-Xmx32G" \
    -sample ${meta.id} \
    -bam_file ${data[0]} \
    -ref_genome ${genome_fasta} \
    -ref_genome_version 38 \
    -ensembl_data_dir ${hmftools_data_dir}/common/ensembl_data/ \
    -threads ${task.cpus} \
    -exp_counts_file ${hmftools_data_dir}/rna/read_151_exp_counts.38.csv \
    -exp_gc_ratios_file ${hmftools_data_dir}/rna/read_100_exp_gc_ratios.38.csv \
    -output_dir isofox/
    """
}











process rseqc {
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/QC/rseqc/", mode: 'copy'


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/rseqc'

    input:
    tuple val(meta), path(data)

    output:
    //path("*.{txt,pdf,r,xls}")
    tuple val(meta), path("*.{txt,pdf,r,xls}")

    when:
    !params.skipQC

    script:
    """
    infer_experiment.py -i ${data[0]} -r ${gene_bed12} > ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc.infer_experiment.txt
    data[0]_stat.py -i ${data[0]} > ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc.data[0]stat.txt
    junction_annotation.py -i ${data[0]} -r ${gene_bed12} -o ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc
    junction_saturation.py -i ${data[0]} -r ${gene_bed12} -o ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc
    inner_distance.py -i ${data[0]} -r ${gene_bed12} -o ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc
    read_distribution.py -i ${data[0]} -r ${gene_bed12} > ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc.read.dist.txt
    read_duplication.py -i ${data[0]} -o ${meta.id}.${meta.npnTumorRNA}_RNA.rseqc.read.dup
    """
}



process qualimapRNAseq {
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/QC/", mode: 'copy'

    input:
    tuple val(meta), path(data)
    
    output:
    path ("${meta.id}.${meta.npnTumorRNA}.qualimapRNA/")

    when:
    !params.skipQC

    script:
    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif qualimap --java-mem-size=20G \
    rnaseq \
    -outdir ${meta.id}.${meta.npnTumorRNA}.qualimapRNA \
    -bam ${bam} -gtf ${gencode_gtf} -pe 
    """
}


process qualimapBAMQC {
    tag "$meta.id"
    publishDir "${meta.id}/QC/RNA/", mode: 'copy'
    cpus 10

    input:
    tuple val(meta), path(data)
    
    output:
    path("${meta.id}.${meta.npnTumorRNA}_RNA.bamqc/")// into qualimapBAMQC
    
    when:
    !params.skipQC && params.qualimap
    
    script:
    //use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    unset DISPLAY
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif qualimap --java-mem-size=5G \
    bamqc \
    -nt ${task.cpus} \
    -outdir ${meta.id}.${meta.npnTumorRNA}_RNA.bamqc \
    -bam ${bam} -gff ${qualimap_ROI} -sd -sdmode 0
    """
}

process rnaseQC {
    tag "$meta.id"
    publishDir "${meta.id}/QC/RNA/RNAseqc", mode: 'copy'
    errorStrategy 'ignore'
    cpus 15

    input:
    tuple val(meta), path(data)

    output:
    path("${meta.npnTumorRNA}_RNA_rnaseqc/*")
    
    when:
    !params.skipQC

    script:
    """
    rnaseqc242 ${gencode_gtf_collapsed} ${bam} \
    --bed ${gencode36_coding_exons} ${meta.npnTumorRNA}_RNA_rnaseqc
    """
}

process multiQC_RNA {
    //publishDir "${meta.id}/", mode: 'copy'
    publishDir "${meta.id}/", mode: 'copy'
    //publishDir "${launchDir}/*/${outputDir}", mode: 'copy'
 
    input:
    //path("*_fastqc.*") from fastqc_results.collect()
    path(inputfiles)
    //path("*.preseq.txt") from preseq_output.collect()
    //path("*.{txt,pdf,r,xls}") from rseqc_out.collect()
    //path("*.{txt,pdf,r,xls}") from rseqc_out.collect()
    //path("*.bamqc/*") from qualimapBAMQC.collect()
    //path("*.qualimapRNA/*") from qualimapRNASEQ.collect()
    //path("*.featureCounts.summary") from featureCounts_output.collect()
    //path("*.RSEM.stat/*") from rsem_stats_out.collect()

    //path("bamQC/*") from bamQCReport.collect()
    output:
    path ("*multiQC.report.html")

    when:
    !params.skipQC

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q ${launchDir}/${meta.id}/QC/DNA/ \
    -n ${date}.RNA_results.multiQC.report.html 
    """
}


process rsem {
    tag "$meta.id"
    cpus 20
    errorStrategy 'ignore'
    
    publishDir "${meta.id}/toolsOutput/genecounts/rsem", mode: 'copy', pattern: "*.results"
    publishDir "${meta.id}/toolsOutput/QC/", mode: 'copy', pattern: "*.${meta.id}.${meta.npnTumorRNA}_RNA.RSEM.stat/*"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/rsem133'

    input:
    tuple val(meta), path(aln)

    output:
    path("*.results")
    tuple val(meta), path("*.genes.results"), emit: rsem_tpm_ch
    path ("*.RSEM.stat/*"),emit: rsem_stats_out
    //    tuple val(meta), path("*pcgrInput*"), emit: rsem_for_pcgr
    script:
 
    strandedness=params.rsem_strand ? "--strandedness ${params.rsem_strand}":""
    """
    rsem-calculate-expression \
    --alignments --no-bam-output -p ${task.cpus} \
    --append-names \
    $strandedness \
    --paired-end ${aln[0]} ${index_rsem} ${meta.id}.${meta.npnTumorRNA}_RNA.RSEM
    """
}

process rsem_genecount_TPM {
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/", mode: 'copy', pattern: "*.2col.txt"
    publishDir "${meta.id}/toolsOutput/genecounts/rsem/", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(meta), path(genecounts)               // from rsem_tpm_ch

    output:
    path("*.tpm.2col.txt"),                     emit: rnaExp2col
    tuple val(meta), path("*.pcgrInput.txt"), emit: rna_for_pcgr
    script:
    """
    cut -f1,6 ${genecounts} > ${meta.npnTumorRNA}_RNA.RSEM.tpm.2col.txt
    sed -i "1 s/TPM/${meta.npnTumorRNA}_RNA/" ${meta.npnTumorRNA}_RNA.RSEM.tpm.2col.txt

    cut -f1,6 ${genecounts} \
    | sed "1 s/gene_id/TargetID/" \
    | awk -F '\t' '{split( \$1, a, ".");print a[1] "\t" \$2 }' > ${meta.npnTumorRNA}_RNA.RSEM.pcgrInput.txt
    """
}


////////////////// ALTERNATIVE SPLICING ///////////////////////////////

process trinitySplicing { 
    errorStrategy 'ignore'
    tag "$meta.id"

    cpus 20
    publishDir "${meta.id}/toolsOutput/trinitySplicing", mode: 'copy'
       // publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.INHOUSE.*"
    input:
    tuple val(meta), path(reads), path(aln), path(starout_sj),path(starout_chimeric) 

    output:
    path("${meta.id}.${meta.npnTumorRNA}_RNA_trinity_splicing.*"), emit: results
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA_trinity_splicing.cancer.introns"), emit: inhouse_list

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/ctat_splicing_v002.sif /usr/local/src/CTAT-SPLICING/STAR_to_cancer_introns.py \
    --SJ_tab_file ${starout_sj} \
    --chimJ_file ${starout_chimeric} \
    --bam_file ${aln[0]} \
    --output_prefix ${meta.id}.${meta.npnTumorRNA}_RNA_trinity_splicing \
    --sample_name ${meta.id}.${meta.npnTumorRNA}_RNA \
    --min_total_reads 10 \
    --vis \
    --ctat_genome_lib ${genome_lib_starfusion}
    """
}

process splicing_inhouse_list {
    errorStrategy 'ignore'
    tag "$meta.id"

    publishDir "${meta.id}/TUMORBOARDFILES/RNA/", mode: 'copy', pattern: "*.INHOUSE.txt"
    publishDir "${meta.id}/toolsOutput/trinitySplicing", mode: 'copy'
    
    input:
    tuple val(meta), path(trinity_splicing)// from inhouse_list_splicing_ch

    output:
    path("*.INHOUSE.txt")

    script:
    """
    cat ${trinity_splicing} | grep -w -f ${inhouse_splicing_genelist} > ${meta.id}.${meta.npnTumorRNA}_RNA.trinity_splicing.INHOUSE.txt
    """
}

///////////////////////////////////////////////////////////////////////////

///////////////////// FUSIONS ///////////////////////////////////

process arriba {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genefusions/arriba", mode: 'copy'

    input:
    tuple val(meta), path(aln)  // from arriba_input_bam_ch //.join(arriba_input_index_ch)
   
    output:
    path("*.{txt,tsv,pdf}"), emit: arriba_out
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.txt"), emit: all_fusions

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/arriba240.sif /arriba_v2.4.0/arriba \
    -x ${aln[0]} -g ${gencode_gtf} -a ${genome_fasta} -b ${arriba_blacklist} -o ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.txt -O ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.discardedfusions.tsv
    
    singularity run -B ${s_bind} ${simgpath}/arriba240.sif /arriba_v2.4.0/draw_fusions.R \
    --fusions=${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.txt \
    --alignments=${aln[0]} \
    --annotation=${gencode_gtf} \
    --cytobands=${arriba_cytoband} \
    --proteinDomains=${arriba_protein_gff} \
    --output=${meta.id}.${meta.npnTumorRNA}_RNA.arribafusionsRplots.pdf

    cut -f1-2 ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.txt | sed "s/\t/--/g" > ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.1col.finspect.txt

    cat ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.txt| grep -w -f ${inhouse_genelist} > ${meta.id}.${meta.npnTumorRNA}_RNA.arriba.fusions.INHOUSEFUSION_V2.txt
    """
}

   

process starfusion {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genefusions/StarFusion/", mode: 'copy'


    input:
    tuple val(meta), path(junctions)

    output:
    path("*.{tsv,txt}"), emit: results
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.StarFusion.INHOUSEFUSION_V2.txt"), emit: fusions 
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA_abridged.tsv"), emit: all_fusions
    
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/ctat_starfusion_v1_11_1.sif STAR-Fusion \
    --genome_lib_dir ${genome_lib_starfusion} \
    -J ${junctions} \
    --examine_coding_effect \
    --output_dir .
    mv star-fusion.fusion_predictions.tsv ${meta.id}.${meta.npnTumorRNA}_RNA_star-fusion.tsv
    mv star-fusion.fusion_predictions.abridged.tsv ${meta.id}.${meta.npnTumorRNA}_RNA_abridged.tsv
    mv star-fusion.fusion_predictions.abridged.coding_effect.tsv ${meta.id}.${meta.npnTumorRNA}_RNA_abridged.coding_effect.tsv

    cat ${meta.id}.${meta.npnTumorRNA}_RNA_abridged.tsv | grep -w -f /data/shared/genomes/databases/genelists/tumortarget/tumortarget.inhouse.v2.127genes.txt > ${meta.id}.${meta.npnTumorRNA}_RNA.StarFusion.INHOUSEFUSION_V2.txt
    """
}

process kallisto_pizzly {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genefusions/pizzly/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/pizzly'
    
    maxForks 5
    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.pizzlyFusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(meta), path("${meta.id}.${meta.npnTumorRNA}_RNA.pizzlyFusions.txt"), emit: all_fusions
    path("${meta.id}_kallisto/")
    path("${meta.id}.${meta.npnTumorRNA}_RNA.fusions.fasta")
    path("${meta.id}.${meta.npnTumorRNA}_RNA.json")
   
    shell:
    '''
    singularity run -B !{s_bind} !{simgpath}/kallisto-0.46.2.sif kallisto quant \
    -i !{kallisto_index} \
    --fusion \
    -o !{meta.id}_kallisto/ \
    !{reads[0]} !{reads[1]}

    pizzly -k 31 \
    --gtf !{ensemble102_gtf} \
    --cache index.cache.txt \
    --align-score 2 \
    --insert-size 400 \
    --fasta !{ensemble102_transcript_fasta} \
    --output !{meta.id}.!{meta.npnTumorRNA}_RNA !{meta.id}_kallisto/fusion.txt
    
    pizzly_flatten_json.py !{meta.id}.!{meta.npnTumorRNA}_RNA.json > !{meta.id}.!{meta.npnTumorRNA}_RNA_pizzlyFusionsRAW.txt
    
    awk 'BEGIN{FS=OFS="\t"} ($5>4||$6>4)' !{meta.id}.!{meta.npnTumorRNA}_RNA_pizzlyFusionsRAW.txt > !{meta.id}.!{meta.npnTumorRNA}_RNA.pizzlyFusions.txt

    cat !{meta.id}.!{meta.npnTumorRNA}_RNA.pizzlyFusions.txt | grep -w -f !{inhouse_genelist} > !{meta.id}.!{meta.npnTumorRNA}_RNA.pizzlyFusions.INHOUSEFUSION_V2.txt
    '''
}

process jaffa_conda {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genefusions/jaffa/", mode: 'copy'

    conda '/data/shared/programmer/miniconda3/envs/bpipe'

    cpus 12
    maxForks 5

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.jaffaFusions.INHOUSEFUSION_V2.csv"), emit: fusions
    tuple val(meta), path("${meta.id}.jaffaFusions.csv"), emit: all_fusions
    path("${meta.id}.jaffaFusions.fasta")
    path("${meta.id}.jaffaFusions.csv")
    script:
    """
    bpipe run \
    -n ${task.cpus} \
    /data/shared/programmer/JAFFA-version-2.3/JAFFA_direct.groovy \
    ${reads[0]} ${reads[1]}
    mv jaffa_results.csv ${meta.id}.jaffaFusions.csv
    mv jaffa_results.fasta ${meta.id}.jaffaFusions.fasta 

    cat ${meta.id}.jaffaFusions.csv | grep -w -f ${inhouse_genelist} > ${meta.id}.jaffaFusions.INHOUSEFUSION_V2.csv
    """
}

process fusioncatcher {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genefusions/", mode: 'copy'
    cpus 30
    maxForks 5
    input:
    tuple val(meta), path(reads)

    output:
    path("fusioncatcher/")
    tuple val(meta), path("fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.INHOUSEFUSION_V2.txt"), emit: fusions
    tuple val(meta), path("fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.txt"), emit: all_fusions
    script:
    """
    singularity run -B ${s_bind} ${simgpath}/fusioncatcher-1.33.sif /opt/fusioncatcher/v1.33/bin/fusioncatcher.py \
    -d ${fusioncatcher_db} \
    -i ${reads[0]},${reads[1]} \
    --skip-blat \
    --skip-star \
    -p ${task.cpus} \
    -o fusioncatcher

    mv fusioncatcher/final-list_candidate-fusion-genes.txt fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.txt

    mv fusioncatcher/final-list_candidate-fusion-genes.vcf fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.vcf

    mv fusioncatcher/summary_candidate_fusions.txt fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.candidateFusionSummary.txt

    cat fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.txt | grep -w -f ${inhouse_genelist} > fusioncatcher/${meta.id}.${meta.npnTumorRNA}_RNA.fusioncatcher.fusions.INHOUSEFUSION_V2.txt
    """
}

process fusionreport_inhouse {
    errorStrategy 'ignore'
    tag "$meta.id"

    publishDir "${meta.id}/toolsOutput/genefusions/fusionReportINHOUSEFUSION_V2/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/", mode: 'copy', pattern: "*.html"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/fusionreport'

    input:
    tuple val(meta), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions), path(pizzly_fusions)//,path(jaffa_fusions)
    output:
    path("fusionreport_output/")
    path("${meta.id}.FusionReport.INHOUSEFUSION_V2.html")
    script:
    """
    fusion_report run \
    ${meta.id} \
    fusionreport_output \
    ${fusionreport_db} \
    --arriba ${arriba_fusions} \
    --starfusion ${starfusion_fusions} \
    --fusioncatcher ${fusioncatcher_fusions} \
    --pizzly ${pizzly_fusions}

    mv fusionreport_output/index.html fusionreport_output/${meta.id}.FusionReport.INHOUSEFUSION_V2.html
    cp fusionreport_output/${meta.id}.FusionReport.INHOUSEFUSION_V2.html .
    """
    //
    //    --jaffa ${jaffa_fusions}
}

process fusionreport_full {
    errorStrategy 'ignore'
    tag "$meta.id"

    publishDir "${meta.id}/toolsOutput/genefusions/fusionReportALL/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/", mode: 'copy', pattern: "*.html"

   conda '/lnx01_data3/shared/programmer/miniconda3/envs/fusionreport'

    input:
    tuple val(meta), path(arriba_fusions), path(starfusion_fusions), path(fusioncatcher_fusions), path(pizzly_fusions)//,path(jaffa_fusions)
    output:
    path("fusionreportALL/")
    path("${meta.id}.FusionReport.ALL.html")
    script:
    """
    fusion_report run \
    ${meta.id} \
    fusionreportALL \
    ${fusionreport_db} \
    --arriba ${arriba_fusions} \
    --starfusion ${starfusion_fusions} \
    --fusioncatcher ${fusioncatcher_fusions} \
    --pizzly ${pizzly_fusions}

    mv fusionreportALL/index.html fusionreportALL/${meta.id}.FusionReport.ALL.html
    cp fusionreportALL/${meta.id}.FusionReport.ALL.html .
    """
}






/*
Arriba:
publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.{INHOUSEFUSION_V2.txt,pdf} "

StarFusion:
    //publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.INHOUSEFUSION_V2.txt"








// Discontinued processes:


process featureCounts {
    tag "$meta.id"
    errorStrategy 'ignore'

    publishDir "${meta.id}/toolsOutput/genecounts/featureCount", mode: 'copy'
   // publishDir "/data/shared/projects/tumortarget/expression/featurecounts/", mode: 'copy'
    publishDir "${meta.id}/QC/", mode: 'copy', pattern: "*.summary"
    cpus 10

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/subread206'
    input:
    tuple val(meta), path(aln)
    
    output:
    path("*.{featureCounts, featureCounts.summary}")
    
    script:
    """
    featureCounts -p -T ${task.cpus} -s 2 -t exon -a ${gencode_gtf} -g gene_id ${aln[0]} -o ${meta.id}.${meta.npnTumorRNA}_RNA.featureCounts
    """
}

process htseq_count {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/genecounts/htseq_count", mode: 'copy'
    //    publishDir "/data/shared/projects/tumortarget/expression/htseq/", mode: 'copy'
    //publishDir "/data/shared/genomes/hg19/databases/rna_seq/inhouse/genecounts/htseq_count/", mode: 'copy'
    cpus 10

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/htseq'
    input:
    tuple val(meta), path(aln)
    
    output:
    path("${meta.id}.${meta.npnTumorRNA}_RNA.htseq_count.txt") 
    
    script:
    """
    htseq-count ${aln[0]} ${gencode_gtf} -s reverse > ${meta.id}.${meta.npnTumorRNA}_RNA.htseq_count.txt
    """
}








    */
    




