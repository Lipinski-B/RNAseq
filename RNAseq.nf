#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Lipinski-B/RNAseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Lipinski-B/RNAseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def nextflowMessage() {
    log.info "N E X T F L O W  ~  DSL 2  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOG INFO
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
log.info "--------------------------------------------------------------------------------------"
log.info ""
log.info "          Pipeline RNAseq for the Differential Expression analysis."
log.info ""
log.info ""


if (params.help) {
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info "  USAGE : nextflow run Lipinski-B/RNAseq --input /data/ --output /output/ "
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info ""
    log.info "nextflow run Lipinski-B/RNAseq [-r vX.X -profile docker/singularity] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--input                       FOLDER                      Folder where you can find your data (fasta/fastq files)."
    log.info "--output                      FOLDER                      Folder where you want to find your result."
    log.info ""
    log.info "Optional arguments:"
    log.info "--mapper                      STRING                      [STAR/BWA] Choose the mapper to use between STAR and BWA MEME (Default : BWA MEM)"
    log.info "--R                           STRING                      [on/off] : Chose to use or not the standard R analyses from the pipeline."
    log.info "--metadata                    FILE                        Path where you can find the XLS file to use as metadata for the R analyse. Mandatory is the option --R in on."
    log.info "--thread                      INT                         Number of thread to use."
  
    exit 0
} else {
    log.info "-------------------------------- Nextflow parameters ---------------------------------"
    log.info ""
    log.info "Project              : $workflow.projectDir"
    log.info "Git repository       : $workflow.repository"
    log.info "Release [Commit ID]  : $workflow.revision [$workflow.commitId]"
    log.info "User Name            : $workflow.userName"
    log.info "Run Name             : $workflow.runName"
    log.info "Resume               : $workflow.resume"
    log.info "Script Name          : $workflow.scriptName"
    log.info "Script File          : $workflow.scriptFile"
    log.info "Home Directory       : $workflow.homeDir"
    log.info "Work Directory       : $workflow.workDir"
    log.info "Launch Directory     : $workflow.launchDir"
    log.info "Command line         : $workflow.commandLine"
    log.info "Config Files         : $workflow.configFiles"
    log.info "Config Profile       : $workflow.profile"
    log.info "Container Engine     : $workflow.containerEngine"
    log.info "Container            : $workflow.container"
    log.info "Session ID           : $workflow.sessionId"
    log.info "Script ID            : $workflow.scriptId"
    log.info ""
    log.info "-------------------------------- Workflow parameters ---------------------------------"
    log.info ""
    log.info "date                 : ${params.date}"
    log.info "input                : ${params.input}"
    log.info "output               : ${params.output}"
    log.info ""
    log.info "--------------------------------------------------------------------------------------"
    log.info ""
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE : Run main Lipinski-B/RNAseq analysis pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    RNAseq(Channel.fromFilePairs("${params.input}/*_R{1,2}.fastq.gz", checkIfExists:true))
}




process KALISTO_INDEX {
    tag "INDEX"
    cpus 20
    memory 40.GB
    publishDir "${params.output}/00_INDEX", mode: 'copy'
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
    publishDir "${params.output}/01_KALLISTO", mode: 'copy'
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
    publishDir "${params.output}/01_KALLISTO", mode: 'copy'
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


workflow.onComplete {
	this.nextflowMessage()
	log.info "Completed at: " + workflow.complete
	log.info "Duration    : " + workflow.duration
	log.info "Success     : " + workflow.success
	log.info "Exit status : " + workflow.exitStatus
	log.info "Error report: " + (workflow.errorReport ?: '-')}

workflow.onError {
	this.nextflowMessage()
	log.info "Workflow execution stopped with the following message:"
	log.info "  " + workflow.errorMessage}