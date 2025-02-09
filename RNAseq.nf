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
log.info " "
log.info " "
log.info "          RNAseq pipeline for the Differential Expression analysis."
log.info " "
log.info " "

if (params.help) {
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info "  USAGE : nextflow run Lipinski-B/RNAseq --input /data/ --output /output/ "
    log.info "------------------------------------------------------------------------------------------------------------------------------"
    log.info " "
    log.info "nextflow run Lipinski-B/RNAseq [-r v1.0 -profile docker] [OPTIONS]"
    log.info " "
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--input                       FOLDER                      Folder where you can find your data (fasta/fastq files)."
    log.info "--output                      FOLDER                      Folder where you want to find your result."
    log.info "--protocol                    STRING                      [PE/SE] Choose the protocol to analyse between Paired-end (PE) and Single-end (SE) (Default : PE)"
    log.info "--metadata                    FILE                        File to use as metadata file"
    log.info ""
    log.info "Optional arguments:"
    log.info ""
    log.info "--cpu                         INT                         Number of CPU to uses (Default : 2)"



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
    log.info "protocol             : ${params.protocol}"
    log.info "cpu                  : ${params.cpu}"
    log.info ""
    log.info "--------------------------------------------------------------------------------------"
    log.info " "
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
        FASTQC(fastq)
        KALISTO_INDEX()

        // Exercice 1 : Perform Differential Expression Analysis at transcript level comparing untreated with treated samples.
        // Exercice 2 : Develop a Nextflow pipeline that performs isoform quantification using Kallisto and assembles the
        //              resulting counts into a matrix format. The pipeline should be capable of switching between single-end and
        //              paired-end modes based on Nextflow parameters.
        // Exercice 2 : Implement the DE analysis with Sleuth or DESeq2 as a Nexflow module, using as input a metadata file
        //              defining the condition of each FASTQ files.

        if (params.protocol == "SE") {
            KALISTO_SINGLE_END(fastq, KALISTO_INDEX.out.KALISTO_INDEX_result)
            DESEQ2_ANALYSIS(KALISTO_SINGLE_END.out.EXRACT_INFO_ch.collect())
            MULTIQC(FASTQC.out.FASTQC_result.collect(), KALISTO_SINGLE_END.out.EXRACT_MULTIQC_ch.collect())
    
        } else if (params.protocol == "PE") {
            KALISTO_PAIRED_END(fastq, KALISTO_INDEX.out.KALISTO_INDEX_result)
            DESEQ2_ANALYSIS(KALISTO_PAIRED_END.out.EXRACT_INFO_ch.collect())
            MULTIQC(FASTQC.out.FASTQC_result.collect(), KALISTO_PAIRED_END.out.EXRACT_MULTIQC_ch.collect())
    
        // Exercice 1 : Compare Single-end and Paired-end based on the number of Differentially expressed transcripts.
        } else {
            KALISTO_SINGLE_END(fastq, KALISTO_INDEX.out.KALISTO_INDEX_result)
            KALISTO_PAIRED_END(fastq, KALISTO_INDEX.out.KALISTO_INDEX_result)
            //DESEQ2_COMPARE() --> Custom Rmd report in bin folder
        }

        
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




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    AVAILABLE PROCESS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQC {
    publishDir "${params.output}/${params.protocol}/QC", mode: 'copy', pattern: '{*fastqc.zip}'
	cpus "${params.cpu}"
        
    input:
        tuple val(ID), path(fastq)
	
    output:
	    path("*.zip"), emit : FASTQC_result

	shell:
	R1 = fastq.find { it =~ /_R1\.fastq\.gz$/ }
    R2 = fastq.find { it =~ /_R2\.fastq\.gz$/ }

    '''
    if [ !{params.protocol} == "PE" ] ; then
        pairs="!{R1} !{R2}"
    else
        pairs="!{R2}"
    fi

	fastqc -t !{task.cpus} ${pairs}
    '''
}

process KALISTO_INDEX {
    tag "INDEX"
    publishDir "${params.output}/${params.protocol}/QUANTIFICATION", mode: 'copy'
    cpus "${params.cpu}"

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
    publishDir "${params.output}/${params.protocol}/QUANTIFICATION", mode: 'copy'
    cpus "${params.cpu}"

    input:
        tuple val(ID), path(fastq)
        path(index)

    output:
        path("*.{tsv,json}"), emit : EXRACT_INFO_ch
        path("kallisto_output_*"),            emit : EXRACT_MULTIQC_ch

    script:
    R2 = fastq.find { it =~ /_R2\.fastq\.gz$/ }
    """
    kallisto quant -l 550 -s 150 -b 10 -i ${index} -o ${ID}/ --single -t ${task.cpus} ${R2} > kallisto_output_${ID}.log 2>&1
    mv ${ID}/abundance.tsv ${ID}_SE_abundance.tsv
    mv ${ID}/run_info.json ${ID}_SE_run_info.json
    """
}

process KALISTO_PAIRED_END {
    tag "${ID}"
    publishDir "${params.output}/${params.protocol}/QUANTIFICATION", mode: 'copy'
    cpus "${params.cpu}"

    input:
        tuple val(ID), path(fastq)
        path(index)

    output:
        path("*.{tsv,json}"), emit : EXRACT_INFO_ch
        path("kallisto_output_*"),            emit : EXRACT_MULTIQC_ch

    script:
    """
    kallisto quant -l 550 -s 150 -b 10 -i ${index} -o ${ID}/ -t ${task.cpus} ${fastq} > kallisto_output_${ID}.log 2>&1
    mv ${ID}/abundance.tsv ${ID}_PE_abundance.tsv
    mv ${ID}/run_info.json ${ID}_PE_run_info.json
    """
}

process MULTIQC {
    cpus '1'
    publishDir "${params.output}/${params.protocol}/QC", mode: 'copy'

    input:
        path(file) 
        path(info)

    output:
        path("multiqc_report.html")
        path("multiqc_data")
    
    shell:
    '''
    multiqc .
    '''
}

process DESEQ2_ANALYSIS {
    publishDir "${params.output}/${params.protocol}/", mode: 'copy'
    stageInMode 'copy'

    input:
        path(info)
        
    output:
        path("*.html")

    script:
    """
        
    Rscript ${baseDir}/bin/report.R ${params.protocol} ${params.metadata} "${params.output}/${params.protocol}/QUANTIFICATION" "${baseDir}" "." #"${params.output}/${params.protocol}"
    
    # Rscript /home/bobo/Bureau/Git/script/RNAseq/bin/report.R SE "/home/bobo/Bureau/Git/script/RNAseq/dependencies/metadata.SE.csv" "/home/bobo/Bureau/Git/result/01_KALLISTO" "/home/bobo/Bureau/Git/script/RNAseq" "/home/bobo/Bureau/Git/result" 
    # Rscript /home/bobo/Bureau/Git/script/RNAseq/bin/report.R PE "/home/bobo/Bureau/Git/script/RNAseq/dependencies/metadata.PE.csv" "/home/bobo/Bureau/Git/result/01_KALLISTO" "/home/bobo/Bureau/Git/script/RNAseq" "/home/bobo/Bureau/Git/result" 
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