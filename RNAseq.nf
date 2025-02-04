#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Lipinski-B/RNAseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Lipinski-B/RNAseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


//
// WORKFLOW: Run main Lipinski-B/RNAseq analysis pipeline
//

workflow {
    RNAseq(Channel.fromFilePairs("/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/data/demultiplexed/*_R{1,2}.fastq.gz", checkIfExists:true))
}


//RUN MAIN WORKFLOW
workflow RNAseq {
    take:
        fastq

    main:    
        KALISTO_INDEX()
        KALISTO_SINGLE_END(fastq,KALISTO_INDEX.out.KALISTO_INDEX_result)
        KALISTO_PAIRED_END(fastq,KALISTO_INDEX.out.KALISTO_INDEX_result)
        
    
    //emit:
        //result = Result_all.out
}



process KALISTO_INDEX {
    tag "INDEX"
    cpus 20
    memory 40.GB
    publishDir "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/result/00_INDEX", mode: 'copy'
    container  "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/script/kallisto.0.51.1.simg" : null }"
    containerOptions "--bind /mnt:/mnt"

    output:
        path("Homo_sapiens.GRCh38.index"), emit : KALISTO_INDEX_result

    script:
    """
    # https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

    wget https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    kallisto index --threads ${task.cpus} -i Homo_sapiens.GRCh38.index Homo_sapiens.GRCh38.cdna.all.fa.gz
    
    """
}



process KALISTO_SINGLE_END {
    tag "${ID}"
    cpus 4
    memory 16.GB
    publishDir "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/result/01_KALLISTO", mode: 'copy'
    container  "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/script/kallisto.0.51.1.simg" : null }"
    containerOptions "--bind /mnt:/mnt"

    input:
        tuple val(ID), path(fastq)
        path(index)

    output:
        path("*.{tsv,json}"), emit : EXRACT_INFO_ch

    script:
    R2 = fastq.find { it =~ /_R2\.fastq\.gz$/ }
    """
    kallisto quant -l 550 -s 150 -b 10 -i ${index} -o ${ID}/ --single -t ${task.cpus} ${R2}
    mv ${ID}/abundance.tsv ${ID}_SE_abundance.tsv
    mv ${ID}/run_info.json ${ID}_SE_run_info.json
    """
}



process KALISTO_PAIRED_END {
    tag "${ID}"
    cpus 4
    memory 16.GB
    publishDir "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/result/01_KALLISTO", mode: 'copy'
    container  "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  "/mnt/datagenetique/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/script/kallisto.0.51.1.simg" : null }"
    containerOptions "--bind /mnt:/mnt"

    input:
        tuple val(ID), path(fastq)
        path(index)

    output:
        path("*.{tsv,json}"), emit : EXRACT_INFO_ch

    script:
    """
    kallisto quant -l 550 -s 150 -b 10 -i ${index} -o ${ID}/ -t ${task.cpus} ${fastq}
    mv ${ID}/abundance.tsv ${ID}_PE_abundance.tsv
    mv ${ID}/run_info.json ${ID}_PE_run_info.json
    """
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
