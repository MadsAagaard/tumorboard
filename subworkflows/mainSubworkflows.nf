#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


include {
        markDup;
        align;
        markAdapters;
        fastq_to_ubam;
        inputFiles_symlinks_fastq;
        haplotypecaller;
        mutect2;
        strelka2;
        strelka2_edits;
        manta_somatic;
        cobalt;
        amber;
        sage;
        purple_pass;
        linx;
        hrd_scores_PASS;
        msisensor;
        sequenza_conda;
        sequenza_R_output;
        pcgr_v212_strelka2;
        pcgr_v212_strelka2_manualFilter;
        pcgr_v212_mutect2;
        mosdepth;
        bamtools;
        qualimap;
        collectWGSmetrics;
        multiQC
         } from "../modules/dnaModules.nf" 

include {inputLinks;
        fastp_TRIM;
        align_STAR;
        align_STAR_hmf;
        rseqc;
        qualimapRNAseq;
        qualimapBAMQC;
        rnaseQC;
        //featureCounts;
        //htseq_count;
        rsem;
        rsem_genecount_TPM;
        isofox;
        trinitySplicing;
        splicing_inhouse_list;
        arriba;
        starfusion;
       // kallisto;
        kallisto_pizzly;
        jaffa_conda;
        fusioncatcher;
        fusionreport_inhouse;
        fusionreport_full;
         } from "../modules/rnaModules.nf" 




workflow DNA_PREPROCESS {

    take:
    case_fastq_input_ch     
   
    main:
    inputFiles_symlinks_fastq(case_fastq_input_ch)
    fastq_to_ubam(case_fastq_input_ch)
    markAdapters(fastq_to_ubam.out)
    align(markAdapters.out)
    markDup(align.out)

    emit:
    finalAln=markDup.out.markDup_output //caseID, sampleID, CRAM, CRAI,type
}

workflow DNA_QC {

    take: 
    dnaInputCRAM_single
    main:
    qualimap(dnaInputCRAM_single)
    collectWGSmetrics(dnaInputCRAM_single)
    mosdepth(dnaInputCRAM_single)

   // collectWGSmetrics.out.multiqc.ifEmpty([]).mix(qualimap.out.multiqc.ifEmpty([])).mix(mosdepth.out.multiqc.ifEmpty([]))

    qualimap.out.multiqc.ifEmpty([]).mix(mosdepth.out.multiqc.ifEmpty([]))
    |groupTuple
    |set {multiqcInput}
    multiqcInput.view()
    multiQC(multiqcInput)


    emit:
    qualimap=qualimap.out.multiqc
    collectWGSmetrics=collectWGSmetrics.out.multiqc
    mosdepth=mosdepth.out.multiqc


}
/*
workflow DNA_QC {

    take: 
    dnaInputCRAM_single
    main:
    qualimap(dnaInputCRAM_single)
    collectWGSmetrics(dnaInputCRAM_single)
    bamtools(dnaInputCRAM_single)
    mosdepth(dnaInputCRAM_single)

    collectWGSmetrics.out.multiqc.ifEmpty([]).mix(qualimap.out.multiqc.ifEmpty([])).mix(mosdepth.out.multiqc.ifEmpty([])).mix(bamtools.out.multiqc.ifEmpty([])).collect()
    |set {multiqcInput}
    multiqcInput.view()

    multiQC(multiqcInput)
}
*/



workflow DNA_STANDARD {

    take:
    dnaInputCRAM
    //rna_for_pcgr
    main:
    haplotypecaller(dnaInputCRAM)
    mutect2(dnaInputCRAM)
    strelka2(dnaInputCRAM)
    strelka2_edits(strelka2.out.strelkarenameVCF)
    msisensor(dnaInputCRAM)
    sequenza_conda(dnaInputCRAM)
    sequenza_R_output(sequenza_conda.out)
    pcgr_v212_strelka2(strelka2_edits.out.strelka2_PASS)
    pcgr_v212_strelka2_manualFilter(strelka2_edits.out.strelka2_PASS_TMB_filtered)

    if (params.skipRNA || params.skipExpression) {
        mutect2.out.mutect2_tumorPASS
        |map {meta,vcf,idx -> tuple(meta, [vcf,idx])}
        |set {mutect2_tumorPASS_edit}

        pcgr_v212_mutect2(mutect2_tumorPASS_edit)
    }

    emit:
    mutect2_PASS=mutect2.out.mutect2_PASS                   //meta,vcf,idx
    mutect2_tumorPASS=mutect2.out.mutect2_tumorPASS         //meta,vcf,idx
    mutect2_snpEFF=mutect2.out.mutect2_snpEFF               //meta,vcf,idx
}



workflow DNA_WGS {

    take:
    dnaInputCRAM
    main:
    manta_somatic(dnaInputCRAM)
    amber(dnaInputCRAM)
    cobalt(dnaInputCRAM)
    sage(dnaInputCRAM)

    amber.out.join(cobalt.out).join(manta_somatic.out.mantaSV_pass).join(sage.out.sage_pass)
    //  | map {meta, amber, cobalt, manta, sage -> tuple(meta, [amber, cobalt, manta,sage ])}
    | set {purple_pass_input}

    purple_pass(purple_pass_input)
    linx(purple_pass.out.purpleDir)


    sage.out.sage_pass.join(purple_pass.out.purple_pass_for_hrd).join(manta_somatic.out.mantaSV_pass)
    //  | map {meta, sage, purple, manta -> tuple(meta, [sage, purple, manta])}

    | set {hrd_PASS_input}

    hrd_scores_PASS(hrd_PASS_input)
    emit:
    cna_for_pcgr=purple_pass.out.cna_for_pcgr
}

workflow RNA_PREPROCESS {

    take:
    rnaInputFASTQ

    main:
    inputLinks(rnaInputFASTQ)
    fastp_TRIM(rnaInputFASTQ)
    align_STAR(fastp_TRIM.out.fastp_out_ch) // star_out_bam_ch: meta.caseID, sampleID, bam,bai
    align_STAR_hmf(fastp_TRIM.out.fastp_out_ch)
    emit:
    star_out_cram_ch=align_STAR.out.star_out_cram_ch // meta.caseID, sampleID, bam, bai
    rsem_input_bam=align_STAR.out.rsem_input_bam  // meta.caseID, sampleID, transcriptomeBAM
    arriba_input_bam=align_STAR.out.arriba_input_bam // Above: meta.caseID, sampleID, bam, bai 
   // star_junctions=align_STAR.out.star_junction_out   // meta.caseID, jounctionfile
    chimeric_junctions_out=align_STAR.out.chimeric_junctions_out // meta.caseID, sampleID, junc.file       
    //star_sjtab=align_STAR.out.star_sjtab_out           // meta.caseID, SJ.out.tab
    sj_chimJunction=align_STAR.out.sj_chimJunction
    isofox_input_cram=align_STAR_hmf.out.isofox_input_cram     

}

workflow RNA_EXPRESSION {

    take:
    rnaInputCRAM
    rnaRsemBAM
    isofoxInputCRAM
    main:
    //featureCounts(rnaInputCRAM)
    //htseq_count(rnaInputCRAM)
    rsem(rnaRsemBAM)
    rsem_genecount_TPM(rsem.out.rsem_tpm_ch)
    isofox(isofoxInputCRAM)
    emit:
    //featureCount_exp=featureCounts.out
    //htseqCount_exp=htseq_count.out
    //rsem_exp=rsem_genecount_TPM.out.rnaExp2col
    //rsem_stats=rsem.out.rsem_stats_out
    rna_for_pcgr=rsem_genecount_TPM.out.rna_for_pcgr
    isofoxDir=isofox.out.isofoxDir
}

workflow RNA_SPLICING {

    take:
    trinity_splicing_input

    main:
    trinitySplicing(trinity_splicing_input)
    splicing_inhouse_list(trinitySplicing.out.inhouse_list)


}

workflow RNA_FUSION {

    take:
    rnaInputFASTQ // raw input channel: meta [r1,r2]
    star_arriba_bam     // meta [bam, bai]
    star_chimeric_junctions // meta [junction file] (for starfusion)

    main:
    arriba(star_arriba_bam)

    fusioncatcher(rnaInputFASTQ)     // TESTING
    kallisto_pizzly(rnaInputFASTQ)   // TESTING
    starfusion(star_chimeric_junctions)
   
    arriba.out.fusions
        .join(starfusion.out.fusions)
        .join(fusioncatcher.out.fusions)
        .join(kallisto_pizzly.out.fusions)
    .set{fusionreport_input}
   
    arriba.out.all_fusions
        .join(starfusion.out.all_fusions)
        .join(fusioncatcher.out.all_fusions)
        .join(kallisto_pizzly.out.all_fusions)
    .set{fusionreport_allFusions_input}
   
    fusionreport_inhouse(fusionreport_input)
    fusionreport_full(fusionreport_allFusions_input)
}

workflow RNA_QC {
    take:
    rnaInputCRAM                     

    main:
    rseqc(rnaInputCRAM)
    qualimapRNAseq(rnaInputCRAM)
    qualimapBAMQC(rnaInputCRAM)
    //rnaseQC(rnaInputCRAM)

    emit:
    rseqc_out=rseqc.out
    qualimapRNA_out=qualimapRNAseq.out
    qualimapBAMQC_out=qualimapBAMQC.out
    //rnaseQC_out=rnaseQC.out
}
