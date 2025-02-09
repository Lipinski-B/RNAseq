# RNA-Pipeline 
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/repository/docker/lipinskiboris/rnaseq)

## Description of my results

### Nextflow pipeline
In this git repository, you will find all the code written to build the nextflow pipeline. This pipeline has been developed to perform RNAseq analyses from FASTQ files and to provides answers to the test exercices

I invite you to read and follow the Dependencies, Parameters and Usage sections if you want to test it out.

Two firsts reports are available in the result folder : report_PE.html + report_SE.html. These reports are the typical standard output you will get when you run the pipeline.

### Biologicals results

Also, you will find a last HTML report that I have been able to produce from my computer in the result folder, to give answers to the exercice 1 and 2, and also to this specific task : Compare Single-end and Paired-end based on the number of Differentially expressed transcripts.

This is the main HTML you should read at first : result/comparison_SE_PE.html

Inside this last HTML, you will be able to find :
  * boxplots to compare kallisto results between protocols Single-End (SE) and Paired-End (PE),
  * all the results from the report_PE.html and report_SE.html files,
  * an intersection of the different differential transcripts expression lists produced to count the number of differentially expressed transcripts between the 2 conditions : Treated and Untreated,
  * and a multivariate experimental model of differential expression to obtain the result by processing all the variability at once (Protocol/Contition). I write it as experimental because I think it would need more time to be right and well developed.

I would have liked so much to take more time to make the workflow cleaner, such as creating the nextflow module, improving protocol switching, making the use of the metadata file more permissive, and further developing the DEA report (Gene ontology enrichment, cleaner figure and code, gene regularoty network, improve the modularity, ... ). 

Also, a shiny app would be better to produce a better report, the R code here includes options to vary padj and FC thresholds, the initial idea was to perform a report by modify dynamically these thresholds. 

Finally, the report is a little light on the interpretation of results, I would have preferred to go deeper with the main results.



## Dependencies
The pipeline is functional under Linux distribution and can be used everywhere.

1. Nextflow:
This pipeline is entirely based on the use of [Nextflow](https://www.nextflow.io) (tested with nextflow 24.10.4). It is strongly recommended to read its [documentation](https://www.nextflow.io/docs/latest/getstarted.html) before running the pipeline.

2. Docker: 
Docker containers have also been developed to allow users to run the pipeline without having to install all the necessary dependencies. Installation of the tools [Docker](https://docs.docker.com/engine/install/ubuntu/) are required beforehand. See the last section of "Usage" for more details. 
The container to use and commands for its installation and use can be found [here](https://hub.docker.com/r/lipinskiboris/rnaseq/). You don't need to download it to run the pipeline, unless the actual git configuration doesn't match with your configuration. In this last case scenario, a modificaiton in the nextflow config file could be necessary.
IMPORTANT : Use Docker as sudo to make it work.

3. Additional information:
    * Place all of your FASTQ in the input data directory.
    * Choose the linked metadata file you want to use, they are available here in the dependencies folder.
    * Differential Expression Analysis (DEA) was only perform with DESeq2
    * The pipeline was not configured to be run with slurm (yet), you can run it locally (tested on : 4cpu/16gb)


## Parameters

| Name         | Value         | Description     |
|--------------|---------------|-----------------|
| --input      | Folder / str  | Path to the folder where the FASTQ files to be used for the analysis are located. Make sure you only have the FASTQ files of interest in this folder and nothing else. |
| --output     | Folder / str  | Path to the folder where the results from the pipeline will be stored. |
| --protocol   | [PE/SE] / str | Choose the protocol to analyse between Paired-end (PE) and Single-end (SE) (Default : PE). |
| --metadata   | File / str    | File to use as metadata file. |
| --cpu        | n / int       | Number of CPU to uses (Default : 2) |


## Usage

- Basic pipeline launch :
```bash
nextflow run Lipinski-B/RNAseq --input /data --output /result -profile docker
```

- Pipeline launch for Single-end analysis :

```bash
nextflow run Lipinski-B/RNAseq --input /data --output /result --protocol "SE" -profile docker
```

- In case the pipeline is updated, it is very usefull to download the last release before to re-run it : 

```bash
nextflow pull Lipinski-B/RNAseq
```

## Exemple of a first use

```bash
WORK=/home/ANALYSIS ; cd $WORK
DATA=$WORK/data
RSLT=$WORK/result
nextflow run Lipinski-B/RNAseq --input $DATA --output $RSLT -profile docker
```




## Contributions

  | Name              | Email                       | Description                               |
  |-------------------|-----------------------------|-------------------------------------------|
  | Lipinski Boris    | lipinskiboris@gmail.com     | Developer to contact for support          |
  
