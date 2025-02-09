DESeq_process  <- function(metadata, counts.files, Protocol_lvl="SE", formule){
  
  design <- as.formula(paste("~", formule))
  
  metadata.SE <- metadata[metadata$protocol == Protocol_lvl,]
  
  tx2gene <- read_csv(paste0(script.dir,"/dependencies/tx2gene.gencode.v27.csv"))
  count.txi.SE <- tximport(counts.files[metadata.SE$sample], type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
  
  dds <- DESeqDataSetFromTximport(count.txi.SE, colData = metadata.SE, design = design)
  smallestGroupSize <- 3
  keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
  dds <- dds[keep,]
  dds$condition <- factor(dds$condition, levels = c("Untreated","Treated"))
  dds <- DESeq(dds)
  return(dds)
}
DESeq_result   <- function(dds){
  res <- results(dds)
  res <- res[!is.na(res$padj),]
  res$threshold <- res$padj < 0.05 ; 
  res[res$log2FoldChange > 0 & res$padj < 0.05, "Expression"] <- "Up" 
  res[res$log2FoldChange < 0 & res$padj < 0.05, "Expression"] <- "Down"
  res[res$padj >= 0.05,  "Expression"] <- "Not Sig"
  res <- res[order(res$padj, decreasing = F),]
  return(res)
}
DT.table       <- function(table, filename = "table"){
  DT::datatable(table, class = 'cell-border stripe', rownames = T, 
                extensions = c("Buttons"), 
                selection = "none",
                options = list(dom = "tipfB",
                               select = list(style = 'os', items = 'row'),
                               buttons = list(
                                 list(extend = 'copy',  filename = filename),
                                 list(extend = 'csv',   filename = filename),
                                 list(extend = 'excel', filename = filename),
                                 list(extend = 'pdf',   filename = filename, exportOptions = list(modifier = list(page = 'current'))),
                                 list(extend = 'print', filename = filename),
                                 'colvis'
                               ))) %>%
    DT::formatSignif(colnames(select_if(table, is.numeric))  , digits=3)
}
Top20          <- function(res, dds, n=20) {
  counts.norm <- counts(dds, normalized = TRUE)
  top_sig     <- as.data.frame(res)
  FC          <- top_sig %>% rownames_to_column(var = "symbol") %>% dplyr::arrange(padj) %>% dplyr::pull(log2FoldChange, symbol) %>% head(n=n) #  Order results by padj values
  top_sig     <- data.frame(counts.norm) %>% rownames_to_column(var = "genes") %>% dplyr::filter(genes %in% names(FC))
  FC          <- FC[order(factor(names(FC), levels = top_sig$genes))]
  top_sig      <- top_sig %>% gather(colnames(top_sig)[2:length(colnames(top_sig))], key = "Sample", value = "counts.norm")
  top_sig$Condition <- rep(dds$condition, each = length(names(FC)))
  top_sig$FC        <- rep(FC,length(unique(top_sig$Sample)))
  ggplot(top_sig) +
    geom_point(aes(reorder(x = genes, FC, max), y = counts.norm, color = Condition), cex = 1.5, position=position_jitter(w=0.1,h=0)) +  
    xlab("Transcripts") + ylab("Normalized Counts") + ggtitle(paste0("Top ",n," Significant")) + theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(plot.title = element_text(hjust = 0.5), axis.text=element_text(size=8))
  
}
Volcano_plot   <- function(res, FC=2, xscale = 6, yscale=200, seuil=0.6, label=NULL) {
  data <- as.data.frame(res)
  color = c("#619CFF","grey","#F8766D")
  filter <- data %>% filter(-log10(data$padj)>1.2 & abs(data$log2FoldChange)>FC)
  filter <- data %>% arrange(padj)
  
  ggplot(data) +
    geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = Expression)) + 
    #ggrepel::geom_label_repel(data=filter[1:20,], aes(label=rownames(data), x = log2FoldChange, y = -log10(padj))) +
    ggtitle("\nVolcano plot :\npadj < 0.05 & log2FoldChange > 1.2\n") + xlab("log2FoldChange") + ylab("-log10 padj") +
    scale_x_continuous(limits = c(-xscale,xscale)) + 
    scale_y_continuous(limits = c(0,yscale)) + 
    scale_color_manual(values = color) + 
    geom_hline(yintercept=1.30103, col="black", linetype="dashed") + 
    geom_vline(xintercept=c(-seuil,seuil), col="black", linetype="dashed") +
    theme_classic() + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1))) 
  
}
MD_plot        <- function(data, xscale = 6, yscale = 5) {
  data <- data[['tr']]
  if(length(unique(data$table$Expression))==2){color = c("grey","#F8766D")}
  else{color = c("#619CFF","grey","#F8766D")}
  ggplot(data$table) +
    geom_point(aes(x = data$AveLogCPM, y = logFC, colour = Expression)) + 
    #ggrepel::geom_label_repel(data=data %>% filter(abs(data$table$logFC)>1.2), aes(label=data$table$genes, x = data$AveLogCPM, y = data$table$logFC)) +
    ggtitle("\nMean-Difference Plot of Pics Data : \nCTRL vs CURCUMIN FDR < 0.05 & logFC > 1.2\n") + xlab("Average log CPM") + ylab("-log10 FC") + 
    scale_x_continuous(limits = c(-1,xscale)) + scale_y_continuous(limits = c(-yscale,yscale)) + scale_color_manual(values = color) + 
    geom_hline(yintercept=c(0.6,-0.6), col="black", linetype="dashed") +
    theme_classic() + theme(plot.title = element_text(size = rel(1.2), hjust = 0.5), axis.title = element_text(size = rel(1))) 
}
MD_plot_DESeq2 <- function(dds){
  par(mfrow=c(2,2),mar=c(2,2,1,1))
  ylim <- c(-2.5,2.5)
  resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
  resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
  resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
  resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
  drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
  plotMA(resGA, ylim=ylim); drawLines()
  plotMA(resLA, ylim=ylim); drawLines()
  plotMA(resG, ylim=ylim); drawLines()
  plotMA(resL, ylim=ylim); drawLines()
}
TopHeatmap     <- function(dds, res, FDR=0.01, FC=0.5){
  annot = as.data.frame(colData(dds)[,c("condition","protocol")])
  data = as.data.frame(res)
  signature <- rownames(data)[data$padj<FDR & abs(data$log2FoldChange)>FC]
  counts.norm <- counts(dds, normalized = TRUE)
  matrix <- counts.norm[signature,]
  pheatmap(matrix, show_rownames = F, scale = "row", main = paste0("Heatmap using normlized read counts\nFDR < ",FDR," & log2FoldChange > ",FC,".\n",nrow(matrix)," targets"), 
           annotation_col = annot, cluster_cols = T)
}
SampleHeatmap  <- function(dds){
  vsd <- vst(dds, blind=FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$protocol, sep=" ")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix, main = "Heatmap sample",
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
}
PCA_plot       <- function(dds){
  vsd <- vst(dds, blind=FALSE)
  pcaData <- plotPCA(vsd, intgroup=c("condition", "protocol"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=condition, shape=protocol)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + theme_bw() + ggtitle("\nPCA\n")
}
Oulier_norm    <- function(dds){
  par(mar=c(8,5,2,2))
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
}
Ontology_Enrichment <- function(condition){
  GO.Up = as.data.frame(topGO(condition,ontology = "BP",sort="down"))
  GO_custom.Up <- data.frame(
    Term = rep(paste(GO.Up$Ont, ":", GO.Up$Term, "\n", GO.Up$N, "genes - FDR =", GO.Up$FDR.Up),2),
    Regulation = c(rep("up",length(GO.Up$Term)),rep("down",length(GO.Up$Term))),
    Count = c(GO.Up$Up,GO.Up$Down), 
    pValueUp = c(GO.Up$FDR.Up, GO.Up$FDR.Up),
    label_ypos = c(GO.Up$Up-3,GO.Up$Up+GO.Up$Down+2)
  )
  plot1 <- ggplot(data=GO_custom.Up, aes(x=reorder(Term, -pValueUp), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7)+ xlab("GO term") + ylab("Number of DEGs up and down inside the up-regulated term.") + ggtitle("Up-regulated pathways")  +
    geom_text(aes(y=label_ypos, label=Count), color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
    theme(aspect.ratio = 3/2, axis.text=element_text(size=17), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  
  GO.Down = as.data.frame(topGO(condition, ontology = "BP" ,sort="up"))
  GO_custom.Down <- data.frame(
    Term = rep(paste(GO.Down$Ont, ":", GO.Down$Term, "\n", GO.Down$N, "genes - FDR =", GO.Down$FDR.Down),2),
    Regulation = c(rep("up",length(GO.Down$Term)),rep("down",length(GO.Down$Term))),
    Count = c(GO.Down$Up,GO.Down$Down), 
    pValueDown = c(GO.Down$FDR.Down, GO.Down$FDR.Down),
    label_ypos = c(GO.Down$Up-3,GO.Down$Up+GO.Down$Down+2)
  )
  plot2 <- ggplot(data=GO_custom.Down, aes(x=reorder(Term, -pValueDown), y=Count, fill=Regulation), color = Regulation) +
    geom_bar(stat="identity", width = 0.7) + xlab("GO term") + ylab("Number of DEGs up and down inside the down-regulated term.") + ggtitle("Down-regulated pathways") +
    geom_text(aes(y=label_ypos, label=Count),  color="black", size=6) + scale_fill_brewer(palette="Paired") + theme_minimal() + 
    theme(aspect.ratio = 3/2, axis.text=element_text(size=15), axis.title=element_text(size=15), plot.title = element_text(size=22)) + coord_flip()
  
  return(list(plot1, plot2))
}
