# RNA-Pipeline 
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/repository/docker/lipinskiboris/rnaseq)

## Description
This pipeline was developed to perform RNAseq analysis from FASTQ files.


## Dependencies
The pipeline is functional under Linux distribution and can be used everywhere.

1. This pipeline is entirely based on the use of [Nextflow](https://www.nextflow.io) (tested with nextflow 24.10.4). It is strongly recommended to read its [documentation](https://www.nextflow.io/docs/latest/getstarted.html) before running the pipeline.

2. Additional information :
  * Place all of your FASTQ in the input data directory.
  * Choose the linked metadata file you want to use, they are available here in the dependencies folder.
  * IMPORTANT : Use Docker as sudo to make it work.
  * Differential Expression Analysis (DEA) was only perform with DESeq2

3. Other: 
Docker containers have also been developed to allow users to run the pipeline without having to install all the necessary dependencies. Installation of the tools [Docker](https://docs.docker.com/engine/install/ubuntu/) are required beforehand. See the last section of "Usage" for more details.


## Results
Here you will find all the results I have been able to produce from my computer. Then, you will find 3 different reports in the git result folder: 

  * report_PE.html & report_SE.html : These reports are the typical results you will get when you run the pipeline. This is a standard output to take the information and give the answers to exercises 1 and 2 of the test.

  * comparison_PE_SE.html : this customized report is designed to provide specific answers to your questions 'Compare Single-end and Paired-end based on the number of Differentially expressed transcript's. Inside, you will be able to find :

      * all the results from the report_PE.html and report_SE.html files

      * an intersection of the different differential gene expression lists produced to count the number of differentially expressed transcripts.

      * and a multivariate experimental model of differential expression to obtain the result by processing all the variability at once. I write it as experimental because I think it would need more time to be right and well developed.

I would have liked so much to take more time to make the workflow cleaner, such as creating the nextflow module, improving protocol switching, making the use of the metadata file more permissive, and further developing the DEA report (Gene ontology enrichment, cleaner figure and code, gene regularoty network, ... ). 

Also, a shiny app would be better to produce a better report, the R code here includes options to vary padj and FC thresholds, the initial idea was to perform a report by modify dynamically these thresholds. 

Finally, the report is a little light on the interpretation of results, I preferred to go straight to the main results.



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
  
