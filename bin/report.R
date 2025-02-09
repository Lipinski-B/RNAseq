#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
protocol <- args[1]  
metadata.file <- args[2] 
result.dir <- args[3] 
script.dir <- args[4]
output <- args[5]


rmarkdown::render(paste0(script.dir,'/bin/report.Rmd'), output_dir='.',
                  params = list(protocol = protocol, metadata.file=metadata.file, result.dir=result.dir, script.dir=script.dir,output=output),
                  #output_file = paste0('report_', protocol, '.html'))
                  output_file = paste0(output,'/report_', protocol, '.html'))

  

# Rscript /home/bobo/Bureau/Git/script/RNAseq/bin/report.R SE "/home/bobo/Bureau/Git/script/RNAseq/dependencies/metadata.SE.csv" "/home/bobo/Bureau/Git/result/01_KALLISTO" "/home/bobo/Bureau/Git/script/RNAseq" "/home/bobo/Bureau/Git/result" 
# Rscript /home/bobo/Bureau/Git/script/RNAseq/bin/report.R PE "/home/bobo/Bureau/Git/script/RNAseq/dependencies/metadata.PE.csv" "/home/bobo/Bureau/Git/result/01_KALLISTO" "/home/bobo/Bureau/Git/script/RNAseq" "/home/bobo/Bureau/Git/result" 