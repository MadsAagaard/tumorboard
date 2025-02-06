#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"

multiqc_config="/data/shared/programmer/configfiles/multiqc_config_tumorBoard.yaml"

//////////////////////////// SWITCHES ///////////////////////////////// 

switch (params.gatk) {

    case 'danak':
    gatk_image="gatk419.sif";
    break;
    case 'new':
    gatk_image="gatk4400.sif";
    break;
    case 'latest':
    gatk_image="gatk4500.sif";
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
    break;
    default:
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/,/fast/:/fast/,/lnx01_data3/:/lnx01_data3/,/lnx01_data4/:/lnx01_data4/";
        simgpath="/data/shared/programmer/simg";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        tmpDIR="/fast/TMP/TMP.${user}/";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/${gatk_image} gatk";
        refFilesDir="/fast/shared/genomes";
    break;
}

switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "/data/shared/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "/data/shared/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "/data/shared/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        break;


    case 'hg38':
        assembly="hg38"
        smncaller_assembly="38"
        pcgr_assembly="grch38"
        // Genome assembly files:
        if (params.hg38v1) {
            genome_fasta = "${refFilesDir}/hg38/GRCh38.primary.fa"
            genome_fasta_fai = "${refFilesDir}/hg38/GRCh38.primary.fa.fai"
            genome_fasta_dict = "${refFilesDir}/hg38/GRCh38.primary.dict"
            genome_version="hg38v1"
            cnvkit_germline_reference_PON="${refFilesDir}/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
            cnvkit_inhouse_cnn_dir="${refFilesDir}/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
            inhouse_SV="${refFilesDir}/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
            genome_fasta = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa"
            genome_fasta_fai = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
            genome_fasta_dict = "${refFilesDir}/hg38/ucsc.hg38.NGS.analysisSet.dict"
            genome_version="hg38v2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "${refFilesDir}/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="hg38v3"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="/data/shared//genomes/hg38/inhouse_DBs/hg38v3_primary/"
        }

        // Program files and resources:
        msisensor_list="/data/shared/genomes/hg38/program_DBs/msisensor/hg38v3_msisensor_scan.txt"
     
        gatk_wgs_pon="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="/data/shared/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

        pcgr_data_dir="/data/shared/genomes/hg38/program_DBs/PCGR/"
        pcgr_VEP="/data/shared/genomes/hg38/program_DBs/PCGRv2/VEP_112_GRCh38_merged/"
        pcgr_data_dir3="/data/shared/genomes/hg38/program_DBs/PCGRv2/20240927/"
        hmftools_data_dir_v534="/data/shared/genomes/hg38/program_DBs/hmftools/v5_34/ref/38"
        hmftools_data_dir_v60="/data/shared/genomes/hg38/program_DBs/hmftools/v6_0/ref/38"
        sequenza_cg50_wig="/data/shared/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"
        virusbreakendDB="/data/shared/genomes/databases/virusbreakenddb_20210401/"
        fakeVirusDir="/data/shared/genomes/hg38/program_DBs/hmftools/fakeVirusDir"
        // Regions & variants:
        qualimap_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="/data/shared/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="/data/shared/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        inhouse127_geneIntervals="/data/shared/genomes/hg38/interval.files/geneIntervals/241022_inhouse127genes.3col.SM.bed"


        callable_regions="/data/shared/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="/data/shared/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        break;
}

// OUTPUT locations:

outputDir="Results"



/*
Default input channel
meta: id, npnNormal,npnTumorDNA,npnTumorRNA,pcgr
data (default): 0:cramNormal, 1:craiNormal, 2:cramTumorDNA, 3:craiTumorDNA, 4:R1_RNA, 5:R2_RNA
*/


process inputFiles_symlinks_cram{
    errorStrategy 'ignore'
    publishDir "${meta.id}/alignments/", mode: 'link', pattern: "*.{ba,cr}*"
    publishDir "${meta.id}/toolsOutput/variantcalls/CramSymlinksDNA/", mode: 'link', pattern: "*.{ba,cr}*"
    //publishDir "${meta.id}/TUMORBOARDFILES/alignments/", mode: 'link', pattern: "*.{ba,cr}*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/alignments/", mode: 'link', pattern: "*.{ba,cr}*"
    input:

    tuple val(meta), path(data)   

    //meta: id, npnNormal,npnTumorDNA,npnTumorRNA,pcgr
    //data (default): 0:cramNormal, 1:craiNormal, 2:cramTumorDNA, 3:craiTumorDNA, 4:R1_RNA, 5:R2_RNA
 
    output:
    tuple val(meta), path("*.symedit.*")
 
    script:
    """
    mv ${data[0]} ${meta.npnNormal}_${params.suffix}.Normal.${genome_version}.symedit.cram
    mv ${data[1]} ${meta.npnNormal}_${params.suffix}.Normal.${genome_version}.symedit.cram.crai
    mv ${data[2]} ${meta.npnTumorDNA}_${params.suffix}.Tumor.${genome_version}.symedit.cram
    mv ${data[3]} ${meta.npnTumorDNA}_${params.suffix}.Tumor.${genome_version}.symedit.cram.crai
    """
}

process inputFiles_symlinks_fastq{
    errorStrategy 'ignore'
    publishDir "${meta.id}/inputSymlinksFASTQ/", mode: 'link', pattern: '*.{fastq,fq}.gz'

    input:

    tuple val(meta), path(data)   
    //meta: id, npnNormal,npnTumorDNA,npnTumorRNA,pcgr
    //data (default): 0:cramNormal, 1:craiNormal, 2:cramTumorDNA, 3:craiTumorDNA, 4:R1_RNA, 5:R2_RNA
    output:
    tuple val(meta), path(data)
    script:
    """
    """
}

process fastq_to_ubam {
    errorStrategy 'ignore'
    tag "$meta.id"
    //publishDir "${outputDir}/unmappedBAM/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${outputDir}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 20
    maxForks 10

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.unmapped.from.fq.bam")
    
    script:
    """
    ${gatk_exec} FastqToSam \
    -F1 ${reads[0]} \
    -F2 ${reads[1]} \
    -SM ${meta.id} \
    -PL illumina \
    -PU KGA_PU \
    -RG KGA_RG \
    --TMP_DIR ${tmpDIR} \
    -O ${meta.id}.unmapped.from.fq.bam
    """
}



process markAdapters {

    input:
    tuple val(meta), path(uBAM)
    
    output:
    tuple val(meta), path("${meta.id}.ubamXT.bam"), path("${meta.id}.markAdapterMetrics.txt")
    
    script:
    """
    ${gatk_exec} MarkIlluminaAdapters \
    -I ${uBAM} \
    -O ${meta.id}.ubamXT.bam \
    --TMP_DIR ${tmpDIR} \
    -M ${meta.id}.markAdapterMetrics.txt
    """
}


process align {
    tag "$meta.id"

    maxForks 6
    errorStrategy 'ignore'
    cpus 60

    input:
    tuple val(meta), path(uBAM), path(metrics)

    output:
    tuple val(meta), path("${meta.id}.${genome_version}.QNsort.BWA.clean.bam")
    
    script:
    """
    ${gatk_exec} SamToFastq \
    -I ${uBAM} \
    -INTER \
    -CLIP_ATTR XT \
    -CLIP_ACT 2 \
    -NON_PF \
    -F /dev/stdout \
    |  singularity run -B ${s_bind} ${simgpath}/bwa0717.sif bwa mem \
    -t ${task.cpus} \
    -p \
    ${genome_fasta} \
    /dev/stdin \
    | ${gatk_exec} MergeBamAlignment \
    -R ${genome_fasta} \
    -UNMAPPED ${uBAM} \
    -ALIGNED /dev/stdin \
    -MAX_GAPS -1 \
    -ORIENTATIONS FR \
    -SO queryname \
    --TMP_DIR ${tmpDIR} \
    -O ${meta.id}.${genome_version}.QNsort.BWA.clean.bam


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        GATK: ${gatk_image}
    END_VERSIONS
    """
}

process markDup {
    errorStrategy 'ignore'
    maxForks 6
    tag "$meta.id"
    publishDir "${meta.id}/alignments/", mode: 'copy', pattern: '*.{ba,cr}*'
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input:
    tuple val(meta), path(aln)
    
    output:
    tuple val(meta),  path("${meta.id}.${genome_version}.BWA.MD.cram"), path("${meta.id}.${genome_version}.BWA.MD*crai"), emit: markDup_output

    script:
    """
    samtools view -h ${aln[0]} \
    | samblaster | sambamba view -t 8 -S -f bam /dev/stdin | sambamba sort -t 8 --tmpdir=${tmpDIR} -o /dev/stdout /dev/stdin \
    |  samtools view \
    -T ${genome_fasta} \
    -C \
    -o ${meta.id}.${genome_version}.BWA.MD.cram -

    samtools index ${meta.id}.${genome_version}.BWA.MD.cram
    """
}

process haplotypecaller {
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 4
    if (params.server=="lnx01") {
    maxForks 2
    }

    publishDir "${meta.id}/toolsOutput/variantcalls/", mode: 'copy', pattern: "*.haplotypecaller.*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.haplotypecaller.*"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}.${meta.npnNormal}.normal.haplotypecaller.vcf.gz"), path("${meta.id}.${meta.npnNormal}.normal.haplotypecaller.vcf.gz.tbi")
    
    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" HaplotypeCaller \
    -I ${data[0]} \
    -R ${genome_fasta} \
    -ERC GVCF \
    -L ${ROI} \
    --smith-waterman FASTEST_AVAILABLE \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --dont-use-soft-clipped-bases \
    -O ${meta.id}.${meta.npnNormal}.normal.g.vcf.gz 
    
    ${gatk_exec} GenotypeGVCFs \
    -R ${genome_fasta} \
    -V ${meta.id}.${meta.npnNormal}.normal.g.vcf.gz  \
    -O ${meta.id}.${meta.npnNormal}.normal.haplotypecaller.vcf.gz \
    -G StandardAnnotation \
    -G AS_StandardAnnotation
    """
}

process mutect2 {
    tag "$meta.id"
    if (params.server=="lnx01") {
    maxForks 2
    }
    
    publishDir "${meta.id}/toolsOutput/variantcalls/", mode: 'copy'
    publishDir "${meta.id}/QC/mutect2_filtering/", mode: 'copy', pattern: '*.{table,stats,tsv}'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: '*.for.VarSeq*'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input:
    tuple val(meta), path(data)
 
    output:

    tuple val(meta), path("${meta.id}.mutect2.for.VarSeq.vcf.gz"), path("${meta.id}.mutect2.for.VarSeq.vcf.gz.tbi"), emit: mutect2_ALL
    tuple val(meta), path("${meta.id}.mutect2.PASSonly.vcf.gz"), path("${meta.id}.mutect2.PASSonly.vcf.gz.tbi"),emit: mutect2_PASS 
    tuple val(meta), path("${meta.id}.mutect2.PASSonly.TUMORonly.vcf.gz"), path("${meta.id}.mutect2.PASSonly.TUMORonly.vcf.gz.tbi"),emit: mutect2_tumorPASS 
    tuple val(meta), path("${meta.id}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz"), path("${meta.id}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz.tbi"), emit: mutect2_PASS_TMB_filtered 
    tuple val(meta), path("${meta.id}.mutect2.PASSonly.snpeff.vcf"), emit: mutect2_snpEFF 
    // for HRD:
    tuple val(meta), path("${meta.id}.mutect2.PASSonly.chr1_22_XY.vcf.gz"), path("${meta.id}.mutect2.PASSonly.chr1_22_XY.vcf.gz.tbi"), emit: mutect2_PASS_reduced

    script:
    """
    ${gatk_exec} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=30" Mutect2 \
    -R ${genome_fasta} \
    -I ${data[2]} \
    -I ${data[0]} \
    -normal ${meta.npnNormal}_${params.suffix} \
    --germline-resource ${mutect_gnomad} \
    --panel-of-normals ${gatk_wgs_pon} \
    -L ${ROI} \
    --dont-use-soft-clipped-bases \
    --native-pair-hmm-threads 30 \
    -pairHMM FASTEST_AVAILABLE \
    --smith-waterman FASTEST_AVAILABLE \
    -O ${meta.id}.mutect2.raw.vcf.gz \
    --f1r2-tar-gz ${meta.id}.within.f1r2.tar.gz
    
    ${gatk_exec} LearnReadOrientationModel \
    -I ${meta.id}.within.f1r2.tar.gz \
    -O ${meta.id}.within.ROmodel.tar.gz
    
    ${gatk_exec} GetPileupSummaries -I ${data[2]} \
    -R ${genome_fasta} \
    -V ${gatk_contamination_ref} \
    -L ${gatk_contamination_ref} \
    -O ${meta.id}.within.getpileupsummaries.table
    
    ${gatk_exec} CalculateContamination \
    -I ${meta.id}.within.getpileupsummaries.table \
    -tumor-segmentation ${meta.id}.segments.table \
    -O ${meta.id}.contamination.table
    
    ${gatk_exec} FilterMutectCalls \
    -V ${meta.id}.mutect2.raw.vcf.gz \
    -R ${genome_fasta} \
    --tumor-segmentation ${meta.id}.segments.table \
    --contamination-table ${meta.id}.contamination.table \
    --min-allele-fraction 0.001 \
    -O ${meta.id}.mutect2.for.VarSeq.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered \
    -O ${meta.id}.mutect2.PASSonly.vcf.gz

    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered \
    -L /data/shared/genomes/hg38/chrom1_22_XY.1col.list \
    -O ${meta.id}.mutect2.PASSonly.chr1_22_XY.vcf.gz
    
    ${gatk_exec} SelectVariants -R ${genome_fasta} \
    -V ${meta.id}.mutect2.for.VarSeq.vcf.gz \
    --exclude-filtered -xl-sn ${meta.npnNormal}_${params.suffix} --exclude-non-variants \
    -O ${meta.id}.mutect2.PASSonly.TUMORonly.vcf.gz

    java -jar /data/shared/programmer/snpEff5.2/snpEff.jar GRCh38.99 ${meta.id}.mutect2.PASSonly.vcf.gz > ${meta.id}.mutect2.PASSonly.snpeff.vcf

    cat ${meta.id}.mutect2.PASSonly.snpeff.vcf | java -jar /data/shared/programmer/snpEff5.2/SnpSift.jar filter \
    "(ANN[0].EFFECT has 'missense_variant'| ANN[0].EFFECT has 'frameshift_variant'| ANN[0].EFFECT has 'stop_gained'| ANN[0].EFFECT has 'conservative_inframe_deletion'|  ANN[0].EFFECT has 'disruptive_inframe_deletion'|ANN[0].EFFECT has 'disruptive_inframe_insertion'|ANN[0].EFFECT has 'conservative_inframe_insertion') & (GEN[${meta.npnTumorDNA}_${params.suffix}].AF >=0.01 & GEN[${meta.npnTumorDNA}_${params.suffix}].DP>25)" > ${meta.id}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf

    bgzip ${meta.id}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    bcftools index -t ${meta.id}.mutect2.PASSonly.snpeff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz
    """
}

process strelka2 {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/variantcalls/strelka2", mode: 'copy'
    cpus 10

    if (params.server=="lnx01"){
            conda '/data/shared/programmer/miniconda3/envs/py310'
    }
    
    input:
    tuple val(meta), path(data)

    output:
    path("*.strelka2.*")
    tuple val(meta), path("${meta.id}.strelka2.merged.vaf.vcf.gz"), emit: strelkarenameVCF

    script:
    def datatype=params.wes ? "": ""
    """
    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif /tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${data[0]} \
    --tumorBam ${data[2]} \
    --referenceFasta  ${genome_fasta} \
    $datatype \
    --runDir strelka

    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif python2 strelka/runWorkflow.py \
    -j ${task.cpus} \
    -m local

    python /data/shared/programmer/VCFpytools/add_vaf_strelka2.py \
    --input strelka/results/variants/somatic.indels.vcf.gz \
    --output ${meta.id}.strelka2.indels.vaf.vcf \
    --variant indel

    python /data/shared/programmer/VCFpytools/add_vaf_strelka2.py \
    --input strelka/results/variants/somatic.snvs.vcf.gz \
    --output ${meta.id}.strelka2.snvs.vaf.vcf \
    --variant snv

    ${gatk_exec} MergeVcfs \
    -I ${meta.id}.strelka2.snvs.vaf.vcf \
    -I ${meta.id}.strelka2.indels.vaf.vcf \
    -O ${meta.id}.strelka2.merged.vaf.vcf.gz
    """
}

process strelka2_edits {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/variantcalls/strelka2", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 

    input: 
    tuple val(meta), path(strelkavcf)

    output:
    
    tuple val(meta), path("${meta.id}.strelka2.for.VarSeq.gz"),path("${meta.id}.strelka2.for.VarSeq.gz.tbi"),emit: strelka2_ALL    
    
    tuple val(meta), path("${meta.id}.strelka2.PASSonly.vcf.gz"),path("${meta.id}.strelka2.PASSonly.vcf.gz.tbi"), emit: strelka2_PASS 
    
    tuple val(meta), path("${meta.id}.strelka2.PASSonly.TUMORonly.vcf.gz"),path("${meta.id}.strelka2.PASSonly.TUMORonly.vcf.gz.tbi"), emit: strelka2_TUMOR_PASS

    tuple val(meta), path("${meta.id}.strelka2.PASSonly.snpeff.vcf"), emit: strelka2_PASS_snpeff
    
    tuple val(meta), path("${meta.id}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz"),path("${meta.id}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz.tbi"), emit: strelka2_PASS_TMB_filtered

    path("*.strelka2.*")
    
    shell:
    '''
    printf "TUMOR !{meta.npnTumorDNA}_!{params.suffix}_TUMOR" >> !{meta.id}.strelka_rename.txt
    printf "\nNORMAL !{meta.npnNormal}_!{params.suffix}_NORMAL" >> !{meta.id}.strelka_rename.txt

    bcftools reheader \
    --samples !{meta.id}.strelka_rename.txt \
    -o !{meta.id}.strelka2.for.VarSeq.gz !{strelkavcf}

    bcftools index -t !{meta.id}.strelka2.for.VarSeq.gz
    
    !{gatk_exec} SelectVariants -R !{genome_fasta} \
    -V !{meta.id}.strelka2.for.VarSeq.gz \
    --exclude-filtered \
    -O !{meta.id}.strelka2.PASSonly.vcf.gz

    !{gatk_exec} SelectVariants -R !{genome_fasta} \
    -V !{meta.id}.strelka2.for.VarSeq.gz \
    --exclude-filtered -xl-sn !{meta.npnNormal}_!{params.suffix}_NORMAL --exclude-non-variants \
    -O !{meta.id}.strelka2.PASSonly.TUMORonly.vcf.gz

    java -jar /data/shared/programmer/snpEff5.2/snpEff.jar GRCh38.99 !{meta.id}.strelka2.PASSonly.vcf.gz > !{meta.id}.strelka2.PASSonly.snpeff.vcf

    cat !{meta.id}.strelka2.PASSonly.snpeff.vcf | java -jar /data/shared/programmer/snpEff5.2/SnpSift.jar filter \
    "(ANN[0].EFFECT has 'missense_variant'| ANN[0].EFFECT has 'frameshift_variant'| ANN[0].EFFECT has 'stop_gained'| ANN[0].EFFECT has 'conservative_inframe_deletion'|  ANN[0].EFFECT has 'disruptive_inframe_deletion'|ANN[0].EFFECT has 'disruptive_inframe_insertion'|ANN[0].EFFECT has 'conservative_inframe_insertion') & (GEN[!{meta.npnTumorDNA}_!{params.suffix}_TUMOR].VAF >=0.05 & GEN[!{meta.npnTumorDNA}_!{params.suffix}_TUMOR].DP>25 & GEN[!{meta.npnNormal}_!{params.suffix}_NORMAL].VAF<0.001)" > !{meta.id}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    
    bgzip !{meta.id}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf
    bcftools index -t !{meta.id}.strelka2.PASSonly.snpEff.snpSift.STDFILTERS_FOR_TMB.v2.vcf.gz
    
    '''
}



process msisensor {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/MSIsensor/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*_msi"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/msisensorPro120'

    if (params.server=="lnx01") {
    maxForks 2
    }
    input: 
    tuple val(meta), path(data)
    output:
    path("*_msi*")
 
    script:
    def datatype=params.wes ? "-e ${ROI}": ""
    """
    msisensor-pro msi \
    -d ${msisensor_list} \
    -n ${data[0]} -t ${data[2]} \
    $datatype \
    -g ${genome_fasta} \
    -o ${meta.id}_msi
    """
}

process sequenza_conda {
    errorStrategy 'ignore'
    tag "$meta.id"
    if (params.server=="lnx01") {
    maxForks 2
    }


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenza30'

    //if (!params.server=="lnx01") {    }
    when:
    !params.skipSequenza

    input:
    tuple val(meta), path(data)
    output:
    tuple val(meta), path("${meta.id}.seqz.final.gz") 
    
    script:
    """
    sequenza-utils bam2seqz \
    -n ${data[0]} -t ${data[2]} \
    --fasta ${genome_fasta} \
    -gc ${sequenza_cg50_wig} \
    -o ${meta.id}.seqz.phase1.gz
    sequenza-utils seqz_binning --seqz ${meta.id}.seqz.phase1.gz \
    -w 50 -o ${meta.id}.seqz.final.gz 
    """
}

process sequenza_R_output {  
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/sequenza/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.*"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenzaEnv'
    
    when:
    !params.skipSequenza
    input:
    tuple val(meta),  path(seqz)

    output:
    //path("sequenza/*")
    path("*.{pdf,txt}")
   // path("*_segments.txt")
    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    t2= sequenza.fit(t1, segment.filter=1e6)
    sequenza.results(sequenza.extract = t1, cp.table = t2, sample.id = "${meta.id}", out.dir = getwd(), CNt.max=1000)
    """
}

/*
process sequenza_R_output {  
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutput/", mode: 'copy', pattern: "sequenza/*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*_{segments,alternative_fit,genome_view}.{txt,pdf}"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sequenzaEnv'
    
    when:
    !params.skipSequenza
    input:
    tuple val(meta),  path(seqz)

    output:
    path("sequenza/*")
    path("*.pdf")
    path("*_segments.txt")
    script:
    """
    #!/usr/bin/env Rscript
    library(sequenza)
    t1 = sequenza.extract("${seqz}",verbose=F)
    t2= sequenza.fit(t1, segment.filter=1e6)
    sequenza.results(sequenza.extract = t1, cp.table = t2, sample.id = "${meta.id}", out.dir = "sequenza", CNt.max=1000)

    cp sequenza/${meta.id}_segments.txt .
    cp sequenza/${meta.id}_alternative_fit.pdf .
    cp sequenza/${meta.id}_genome_view.pdf .
    """
}
*/

process pcgr_v212_strelka2 {
    errorStrategy 'ignore'
    publishDir "${meta.id}/toolsOutput/PCGR212/strelka2/", mode: 'copy', pattern: "*.pcgr.*"
    //publishDir "${meta.id}/${params.outdir}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   
    input:
 //   tuple val(meta),  path(data)
    tuple val(meta),  path(vcf), path(idx)
    output:
    path("*.pcgr.*.{xlsx,tsv,html}")
    
    script:
    
    def datatype=params.wes ? "--assay WES" : "--assay WGS"
    """
    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${meta.id}_strelka2 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display missense_only \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${meta.pcgr} \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures
    """
}

process pcgr_v212_strelka2_manualFilter {
    errorStrategy 'ignore'
    publishDir "${meta.id}/toolsOutput/PCGR212/strelka2_manual", mode: 'copy', pattern: "*.pcgr.*"
    //publishDir "${meta.id}/${params.outdir}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   
    input:

    //   tuple val(meta),  path(data)
    tuple val(meta),  path(vcf), path(idx)

    output:
    path("*.pcgr.*.{xlsx,tsv,html}")
    
    script:
    
    def datatype=params.wes ? "--assay WES": "--assay WGS"
    """
    pcgr \
    --input_vcf ${vcf} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${meta.id}_strelka2_manual \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display missense_only \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    --tumor_site ${meta.pcgr} \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures
    """
}

process pcgr_v212_mutect2 {
    errorStrategy 'ignore'
    publishDir "${meta.id}/toolsOutput/PCGR212/mutect2/", mode: 'copy', pattern: "*.pcgr.*"
    publishDir "${meta.id}/TUMORBOARDFILES/",mode: 'copy', pattern:"*.html"
    //publishDir "${meta.id}/${params.outdir}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.flexdb.html"
   
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/pcgr212'  
   


    input:
    tuple val(meta),  path(data)        
    // channel if default (not RNA skipped): meta, [vcf, idx, rnaExp, cna] 
   //  channel if RNA skipped: meta[vcf,idx,cna]

    output:
    path("*.pcgr.*")
    
    script:

    def rnaexp=!params.skipRNA                   ? "--input_rna_expression ${data[2]}" : ""
    def datatype=params.wes                     ? "--assay WES" : "--assay WGS"
    """

    pcgr \
    --input_vcf ${data[0]} \
    --refdata_dir  ${pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${pcgr_VEP} \
    --genome_assembly ${pcgr_assembly} \
    --sample_id ${meta.id}_pcgr212 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display coding_non_silent \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    $datatype \
    $rnaexp \
    --pcgrr_conda /lnx01_data3/shared/programmer/miniconda3/envs/pcgrr212 \
    --estimate_signatures \
    --tumor_site ${meta.pcgr}
    """
}
    /*
        if (!params.wes && !params.skipRNA) {
        def cnaInput="--input_cna ${data[3]}" 
        }

        if (!params.wes && params.skipRNA) {
        def cnaInput="--input_cna ${data[2]}" 
        }
        if (params.wes) {
        def cnaInput="" 
        }
        */


//////////// WGS only modules:

process manta_somatic {
    errorStrategy 'ignore'
    tag "$meta.id"

    publishDir "${meta.id}/toolsOutput/structuralVariants/manta_somatic/", mode: 'copy'
    //publishDir "${outputDir}/structuralVariants/manta/", mode: 'copy', pattern: "*.{AFanno,filtered}.*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.somaticSV.vcf.*"
    cpus 12

    if (params.server=="lnx01") {
    maxForks 1
    }

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 
    input: 
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.manta.*.{vcf,vcf.gz,gz.tbi}")
    tuple val(meta), path("${meta.id}.manta.somaticSV.vcf.gz"), emit: mantaSV_all
    tuple val(meta), path("${meta.id}.manta.somaticSV.PASSonly.vcf.gz"), emit: mantaSV_pass

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif configManta.py \
    --normalBam ${data[0]} \
    --tumorBam ${data[2]} \
    --referenceFasta ${genome_fasta} \
    --callRegions ${manta_callable_regions} \
    --runDir manta

    singularity run -B ${s_bind} ${simgpath}/manta1.6_strelka2.9.10.sif ./manta/runWorkflow.py -j ${task.cpus}

    mv manta/results/variants/candidateSmallIndels.vcf.gz \
    ${meta.id}.manta.candidateSmallIndels.vcf.gz
    
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    ${meta.id}.manta.candidateSmallIndels.vcf.gz.tbi
    
    mv manta/results/variants/candidateSV.vcf.gz \
    ${meta.id}.manta.candidateSV.vcf.gz
    
    mv manta/results/variants/candidateSV.vcf.gz.tbi \
    ${meta.id}.manta.candidateSV.vcf.gz.tbi

    mv manta/results/variants/diploidSV.vcf.gz \
    ${meta.id}.manta.diploidSV.vcf.gz
    
    mv manta/results/variants/diploidSV.vcf.gz.tbi \
    ${meta.id}.manta.diploidSV.vcf.gz.tbi

    gzip -dc ${meta.id}.manta.diploidSV.vcf.gz > ${meta.id}.manta.diploidSV.vcf

    mv manta/results/variants/somaticSV.vcf.gz \
    ${meta.id}.manta.somaticSV.vcf.gz
    
    mv manta/results/variants/somaticSV.vcf.gz.tbi \
    ${meta.id}.manta.somaticSV.vcf.gz.tbi


    bcftools view \
    -i 'FILTER="PASS" | FILTER="."' \
    ${meta.id}.manta.somaticSV.vcf.gz > ${meta.id}.manta.somaticSV.PASSonly.vcf

    bgzip ${meta.id}.manta.somaticSV.PASSonly.vcf
    bcftools index -t ${meta.id}.manta.somaticSV.PASSonly.vcf.gz


    """
}

process cobalt {
    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$meta.id"

    if (params.server=="lnx01") {
    maxForks 1
    }
    cpus 16
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/hmftools/'

    input: 
    tuple val(meta), path(data)

    output: 
    tuple val(meta), path("cobalt/"), emit: cobaltDir
    
    script:
    """
    cobalt "-Xmx16G" \
    -reference ${meta.npnNormal}_${params.suffix} \
    -reference_bam ${data[0]} \
    -tumor ${meta.npnTumorDNA}_${params.suffix} \
    -tumor_bam ${data[2]} \
    -ref_genome ${genome_fasta} \
    -output_dir cobalt \
    -threads ${task.cpus} \
    -gc_profile ${hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp
    """
}

process amber {
    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$meta.id"

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/hmftools/'

    cpus 12
    
    input: 
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("amber/"), emit: amberDir
    
    script:
    """
    amber "-Xmx16G" \
    -reference ${meta.npnNormal}_${params.suffix} \
    -reference_bam ${data[0]} \
    -tumor ${meta.npnTumorDNA}_${params.suffix} \
    -tumor_bam ${data[2]} \
    -ref_genome ${genome_fasta} \
    -output_dir amber \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -loci ${hmftools_data_dir_v534}/copy_number/GermlineHetPon.38.vcf.gz
    """
}


process sage {
    publishDir "${meta.id}/toolsOutput/hmftools/sage/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$meta.id"

    if (params.server=="lnx01") {
    maxForks 2
    }
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 
    cpus 16
    input: 
    tuple val(meta), path(data)

    output:
    tuple val(meta), path("${meta.id}.sage.somatic.pass.vcf.gz"),emit:sage_pass
    tuple val(meta), path("${meta.id}.sage.somatic.all.vcf.gz"),emit:sage_all
    script:
    """
    java -jar /data/shared/programmer/hmftools/sage_v3.4.4.jar \
    -reference ${meta.npnNormal}_${params.suffix} \
    -reference_bam ${data[0]} \
    -tumor ${meta.npnTumorDNA}_${params.suffix} \
    -tumor_bam ${data[2]} \
    -ref_genome ${genome_fasta} \
    -ref_genome_version 38 \
    -threads ${task.cpus} \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -high_confidence_bed ${hmftools_data_dir_v534}/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz \
    -hotspots ${hmftools_data_dir_v534}/variants/KnownHotspots.somatic.38.vcf.gz \
    -panel_bed ${hmftools_data_dir_v534}/variants/ActionableCodingPanel.38.bed.gz \
    -output_vcf ${meta.id}.sage.somatic.all.vcf.gz


    bcftools view \
    -i 'FILTER="PASS" | FILTER="."' \
    ${meta.id}.sage.somatic.all.vcf.gz > ${meta.id}.sage.somatic.pass.vcf

    bgzip ${meta.id}.sage.somatic.pass.vcf
    bcftools index -t ${meta.id}.sage.somatic.pass.vcf.gz

    """
}

process purple_pass {

    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 12

    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.purity.tsv"
        publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.segment.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.{qc,png}"
    //publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.circos.png"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.driver.catalog.somatic.tsv"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/'

    input:
    tuple val(meta), path(amber), path(cobalt), path(manta_sv), path(sage)

    output:
    tuple val(meta), path("purple/"), emit: purpleDir
    tuple val(meta), path("purple/${meta.id}.purple.cnv.somatic.tsv"), emit: purple_pass_for_hrd
    path("${meta.id}.*")
    tuple val(meta), path("${meta.id}.purple.cnv.somatic.forPCGR.txt"), emit:cna_for_pcgr
    script:
    """
    purple "-Xmx8G" \
    -reference ${meta.npnNormal}_${params.suffix} \
    -tumor ${meta.npnTumorDNA}_${params.suffix} \
    -ref_genome ${genome_fasta} \
    -output_dir purple \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -somatic_sv_vcf ${manta_sv} \
    -somatic_vcf ${sage} \
    -gc_profile ${hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -amber ${amber} \
    -cobalt ${cobalt} \
    -driver_gene_panel ${hmftools_data_dir_v534}/common/DriverGenePanel.38.tsv \
    -somatic_hotspots ${hmftools_data_dir_v534}/variants/KnownHotspots.somatic.38.vcf.gz \
    -circos /lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/bin/circos

    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.cnv.somatic.tsv purple/${meta.id}.purple.cnv.somatic.tsv
    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.qc purple/${meta.id}.purple.qc
    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.purity.tsv purple/${meta.id}.purple.purity.tsv
    mv purple/plot/${meta.npnTumorDNA}_${params.suffix}.circos.png purple/${meta.id}.purple.PASS.circos.png

    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.sv.vcf.gz purple/${meta.id}.purple.sv.vcf.gz
    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.sv.vcf.gz.tbi purple/${meta.id}.purple.sv.vcf.gz.tbi

    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.somatic.vcf.gz purple/${meta.id}.purple.somatic.vcf.gz
    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.somatic.vcf.gz.tbi purple/${meta.id}.purple.somatic.vcf.gz.tbi
    
    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.driver.catalog.somatic.tsv purple/${meta.id}.purple.driver.catalog.somatic.tsv

    mv purple/${meta.npnTumorDNA}_${params.suffix}.purple.segment.tsv purple/${meta.id}.purple.segment.tsv

    cp purple/${meta.id}.* .

    cut -f 1,2,3,15,16 ${meta.id}.purple.cnv.somatic.tsv > ${meta.id}.purple.cnv.somatic.forPCGR.txt
    sed -i "1 s/majorAlleleCopyNumber/nMajor/;s/minorAlleleCopyNumber/nMinor/;s/chromosome/Chromosome/;s/start/Start/;s/end/End/" ${meta.id}.purple.cnv.somatic.forPCGR.txt
    """

}

process linx {
  errorStrategy 'ignore'
    tag "$meta.id"
    cpus 12

    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.{tsv,png}"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/linx20beta/'

    input:
    tuple val(meta), path(purple)

    output:
    tuple val(meta), path("linx/")
    script:
    """
    linx "-Xmx8G" \
    -sample ${meta.id} \
    -output_dir linx \
    -ref_genome_version 38 \
    -purple_dir ${purple} \
    -ensembl_data_dir ${hmftools_data_dir_v534}/common/ensembl_data/ \
    -known_fusion_file ${hmftools_data_dir_v534}/sv/known_fusion_data.38.csv \
    -driver_gene_panel ${hmftools_data_dir_v534}/common/DriverGenePanel.38.tsv
    """
}
//-sample ${meta.npnTumorDNA}_${params.suffix} \

process hrd_scores_PASS {
    errorStrategy 'ignore'
    tag "$meta.id"
  
    publishDir "${meta.id}/toolsOutput/HRD/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/",mode: 'copy', pattern: "*.txt"
    cpus 4
    maxForks 3
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sigrap011/'

    input:
    tuple val(meta),path(sage_snv),path(purple_cnv),path(manta_sv)

    output:
    tuple val(meta), path("${meta.id}.HRD_SCORES.txt")
    path("*.txt")
    script:
    """
    #!/usr/bin/env Rscript
    library(sigrap)
    library(BSgenome.Hsapiens.UCSC.hg38)
 
    chord=sigrap::chord_run(vcf.snv="${sage_snv}",vcf.sv="${manta_sv}", sv.caller="manta", sample.name="${meta.id}", ref.genome="hg38")
    hrdetect=sigrap::hrdetect_run(snvindel_vcf="${sage_snv}",sv_vcf="${manta_sv}", nm="${meta.id}", cnv_tsv="${purple_cnv}", genome="hg38")
   
    c1=as.data.frame(chord[2])
    c2=c1[,1:6]
    names(c2)=c("sample","CHORD_p_hrd","CHORD_hr_status","CHORD_hrdtype","CHORD_pBRCA1", "CHORD_pBRCA2")

    d1=data.frame(hrdetect[1], hrdetect[2])
    names(d1)=c("sample","HRDETECT_p_hrd")

    m1=merge(c2,d1,by="sample")
    if (m1[["HRDETECT_p_hrd"]]>0.7) m1[["HRDETECT_verdict"]]="HR_DEFICIENT" else m1[["HRDETECT_verdict"]]="HR_proficient"
    write.table(m1,file="${meta.id}.HRD_SCORES.txt",sep="\t",quote=F,row.names=F)

    """
}

process virus_breakend {
  errorStrategy 'ignore'
    tag "$meta.id"
    cpus 12

    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.{tsv,png}"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/gridss.2.13.2/'

    input:
    tuple val(meta), path(data)     // meta, data[0:cramNormal,1:craiNormal,2:cramTumor,3:craiTumor]

    output:
    tuple val(meta), path("${meta.id}.virusbreakend.vcf"),emit: virusbreakendVCF
    script:
    """
    virusbreakend \
    -r ${genome_fasta} \
    --db ${virusbreakendDB} \
    -o ${meta.id}.virusbreakend.vcf \
    ${data[2]}

    """
}


process virus_interpreter {
  errorStrategy 'ignore'
    tag "$meta.id"
    cpus 12

    publishDir "${meta.id}/toolsOutput/hmftools/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.{tsv,png}"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/gridss.2.13.2/'

    input:
    tuple val(meta), path(data)     // meta, data[0:cramNormal,1:craiNormal,2:cramTumor,3:craiTumor]

    output:
    tuple val(meta), path("${meta.id}.virusbreakend.vcf"),emit: virusbreakendVCF
    script:
    """
    virusbreakend \
    -r ${genome_fasta} \
    --db ${virusbreakendDB} \
    -o ${meta.id}.virusbreakend.vcf \
    ${data[2]}

    """
}

process cuppa {
  errorStrategy 'ignore'
    tag "$meta.id"
    cpus 12

    publishDir "${meta.id}/toolsOutput/hmftools/cuppa/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.{tsv,png}"
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/gridss.2.13.2/'

    input:
    tuple val(meta), path(data)     // meta, data[0:cramNormal,1:craiNormal,2:cramTumor,3:craiTumor]

    output:
    tuple val(meta), path("cuppa/"),    emit: cuppaDir
    script:
    """
    mkdir cuppa/
    java -cp /data/mmaj/downloads/hmftools/cuppa_v2.3.1-rc.3.jar com.hartwig.hmftools.cup.prep.CuppaDataPrep \
    -sample ${meta.id} \
    -ref_genome_version V38 \
    -output_dir cuppa/ \
    -ref_alt_sj_sites /data/shared/genomes/hg38/program_DBs/hmftools/hmf_pipeline_resources.38_v6.0/misc/cuppa/alt_sj.selected_loci.38.tsv.gz \
    -isofox_dir isofox/ \
    -purple_dir purple/purple/ \
    -linx_dir linx/ \
    -virus_dir ${fakeVirusDir} \
    -categories ALL
    """
}
      



////////////// QC modules:

process mosdepth {
    publishDir "${meta.id}/QC/DNA/mosdepth/", mode: 'copy'
    errorStrategy 'ignore'
    tag "$meta.id"

    if (params.server=="lnx01") {
    maxForks 2
    }
    conda '/lnx01_data3/shared/programmer/miniconda3/envs/mosdepth/' 
    cpus 8

    input: 
    tuple val(meta), path(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.npn}_${meta.sampletype}.*")
    tuple val(meta), path("*.global.dist.txt"),emit:multiqc
    script:
    """
    mosdepth \
    -t ${task.cpus} \
    -f ${genome_fasta} \
    --by ${ROI} \
    ${meta.npn}_${meta.sampletype} \
    ${data[0]}

    """


}

process bamtools {
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/QC/DNA/bamtools/", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sambamvcftools/' 
    
    input:
    tuple val(meta),  path(data)
    output:
    tuple val(meta), path("${meta.id}_${meta.npn}_${meta.sampletype}.bamtools.txt"), emit: multiqc

    script:
    """
    samtools view -hb -T ${genome_fasta} ${data[0]} | \
    bamtools stats \
    -in /dev/stdin \
    -insert > ${meta.id}_${meta.npn}_${meta.sampletype}.bamtools.txt
  
    """
  /*cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bamtools: bamtools --version| grep "bamtools" | sed 's/bamtools //g'
    END_VERSIONS
    */
}

process qualimap {
    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 10
    maxForks 6
    publishDir "${meta.id}/QC/DNA/", mode: 'copy'

    conda '/lnx01_data3/shared/programmer/miniconda3/envs/qualimapSamtools/' 
    input:
    tuple val(meta), path(data)

    output:
    tuple val(meta), path ("${meta.npn}_${meta.sampletype}_qualimap"), emit: multiqc
    tuple val(meta), path ("${meta.npn}_${meta.sampletype}_qualimap/qualimapReport.html")
    //path ("versions.yml"), emit: versions

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    samtools view -hb -T ${genome_fasta} ${data[0]} > ${meta.npn}.bam

    qualimap --java-mem-size=40G bamqc \
    --collect-overlap-pairs \
    -nt ${task.cpus} \
    -outdir ${meta.npn}_${meta.sampletype}_qualimap \
    -bam ${meta.npn}.bam \
    $use_bed \
    -sd \
    -sdmode 0
    """
}

/*
process qualimap {
    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 10
    maxForks 8
    publishDir "${meta.id}/QC/DNA/", mode: 'copy'

    input:
    tuple val(meta), path(data)

    output:
    path ("qualimap/"), emit: multiqc
    //path ("versions.yml"), emit: versions

    script:
    use_bed = qualimap_ROI ? "-gff ${qualimap_ROI}" : ''
    """
    unset DISPLAY
    samtools view -hb -T ${genome_fasta} ${data[0]} |
    singularity run -B ${s_bind} ${simgpath}/qualimap.sif \
    qualimap --java-mem-size=5G bamqc \
    --collect-overlap-pairs \
    -nt ${task.cpus} \
    -outdir qualimap \
    -bam  /dev/stdin \
    -sd \
    -sdmode 0
    """
}
*/




process collectWGSmetrics {

    errorStrategy 'ignore'
    tag "$meta.id"
    cpus 5
    publishDir "${meta.id}/QC/DNA/picardWGSmetrics/", mode: 'copy'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.id}_${meta.npn}_${meta.sampletype}.picardWGSmetrics.txt"), emit: multiqc

    when:
    params.picard

    script:
    """
    ${gatk_exec} CollectWgsMetrics \
    -I ${data[0]} \
    -O ${meta.id}_${meta.npn}_${meta.sampletype}.picardWGSmetrics.txt \
    -R ${genome_fasta}
    """
}


process multiQC {
    
    errorStrategy 'ignore'
    publishDir "${meta.id}/QC/", mode: 'copy'

    input:
    //path(inputfiles)
    tuple val(meta),  path(data)  
    //   path("_fastqc.*").collect().ifEmpty([])
    // path("${meta.id}.samtools.sample.stats.txt").collect().ifEmpty([])
    // path("bamQC/*").collect().ifEmpty([]) 
    //path("${meta.id}.picardWGSmetrics.txt").collect().ifEmpty([]) 

    output:
    path ("MultiQC*.html")

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q  ${launchDir}/*/QC/DNA/ \
    -n MultiQC.DNA.html
    """
}

process multiQC_ALL {
    
    errorStrategy 'ignore'
    publishDir "${meta.id}/", mode: 'copy'

    input:
    tuple val(meta),  path(data)  
    //   path("_fastqc.*").collect().ifEmpty([])
    // path("${meta.id}.samtools.sample.stats.txt").collect().ifEmpty([])
    // path("bamQC/*").collect().ifEmpty([]) 
    //path("${meta.id}.picardWGSmetrics.txt").collect().ifEmpty([]) 

    output:
    path ("${meta.id}.MultiQC.ALL.html")

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/multiqc.sif \
    -c ${multiqc_config} \
    -f -q  ${launchDir}/${meta.id}/ \
    -n ${meta.id}.MultiQC.ALL.html
    """
}

/*
qualimap end versions

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Qualimap: v.2.2.1' )
    END_VERSIONS
*/
