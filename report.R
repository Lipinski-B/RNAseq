#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
protocol <- args[1]  
metadata.file <- args[2] 

rmarkdown::render('/home/bobo/Bureau/Git/script/RNAseq/report.Rmd',
                  params = list(protocol = protocol, metadata.file=metadata.file),
                  output_file = paste0('/home/bobo/Bureau/Git/script/RNAseq/report_', protocol, '.html'))

# Rscript report.R SE "/home/bobo/Bureau/Git/script/RNAseq/metadata.SE.csv"