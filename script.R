library(jsonlite)
library(dplyr)
library(ggplot2)
library(tidyverse)
library("DESeq2")

logs.path <- "/home/bobo/Bureau/Git/result/01_KALLISTO/"
logs.files <- list.files(logs.path, pattern = "_*E_run_info.json$", full.names = TRUE)
logs.result <- data.frame()


for (log in logs.files) {
  log.json <- t(as.data.frame(jsonlite::fromJSON(log)))
  log.json <- as.data.frame(log.json[c("n_targets", "n_processed", "n_pseudoaligned", "n_unique", "p_pseudoaligned", "p_unique"), ])
  colnames(log.json) <- gsub("_run_info.json", "", basename(log))
  if (nrow(logs.result) == 0) {
    logs.result <- log.json
  } else {
    logs.result <- cbind(logs.result, log.json)
  }
}

logs.result <- as.data.frame(t(logs.result))
logs.result$Sample <- rownames(logs.result)
logs.result$Condition <- "Traited"
logs.result[grep("Untreated", logs.result$Sample, ignore.case = TRUE), "Condition" ] <- "Untreated"
logs.result$Protocol <- "PE"
logs.result[grep("SE", logs.result$Sample, ignore.case = TRUE), "Protocol" ] <- "SE"

logs.pivot <- logs.result %>% pivot_longer(cols = -c(Sample,Condition,Protocol), names_to = "Feature", values_to = "Value")
logs.pivot$Value <- as.numeric(logs.pivot$Value)

for (col in c("n_targets", "n_processed", "n_pseudoaligned", "n_unique", "p_pseudoaligned", "p_unique")) {
  logs.subset <- logs.pivot %>% filter(Feature == col)
  print(ggplot(logs.subset) +
          aes(x = Protocol, y = Value) + 
          geom_boxplot() +
          scale_fill_hue(direction = 1) +
          theme_minimal() +
          facet_wrap(vars(Feature), scales = "free"))
}


moyennes <- logs.subset %>% group_by(Protocol) %>% summarise(moyenne = mean(Value, na.rm = TRUE))
ggplot(logs.subset) + aes(x = Protocol, y = Value, fill = Condition) +
  geom_boxplot() +
  geom_hline(data = moyennes, aes(yintercept = moyenne, color = Protocol), linetype = "dashed", size = 1) +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  facet_wrap(vars(Feature), scales = "free") +
  labs(color = "Moyenne par protocole")












counts.path <- "/home/bobo/Bureau/Git/result/01_KALLISTO/"
counts.files <- list.files(counts.path, pattern = "_*E_abundance.tsv", full.names = TRUE)
counts.result <- data.frame()


read_kallisto_counts <- function(file) {
  df <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  counts <- df$est_counts
  names(counts) <- df$target_id
  return(counts)
}

count_list <- lapply(counts.files, read_kallisto_counts)
count_matrix <- do.call(cbind, count_list)
count_matrix <- as.data.frame(lapply(count_matrix, function(x) as.numeric(as.character(x))))
count_matrix <- round(count_matrix)
colnames(count_matrix) <- c("Treated1_PE"  , "Treated1_SE"  ,"Treated2_PE" ,  "Treated2_SE"  ,"Treated3_PE"  , "Treated3_SE"  ,"Treated4_PE"   ,"Treated4_SE"  ,"Untreated1_PE" ,"Untreated1_SE","Untreated2_PE" ,"Untreated2_SE","Untreated3_PE" ,"Untreated3_SE","Untreated4_PE" ,"Untreated4_SE")

#count <- count[count$tpm!=0,]

metadata <- data.frame(condition = factor(c(rep("Treated",8), rep("Untreated",8))))
y <- DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = ~condition)
y <- DESeq(y)#, minReplicatesForReplace=Inf)
y$samples$Condition <- y$condition

tr <- result <- list()
tr$table <- results(y, contrast=c("condition",unique(as.vector(y$samples$Condition))[2],unique(as.vector(y$samples$Condition))[1]), pAdjustMethod = "BH", alpha=0.05)#, cooksCutoff=FALSE, independentFiltering=FALSE)
tr$table <- tr$table[complete.cases(tr$table),]
tr$table <- tr$table[order(tr$table$padj),]
tr$table$threshold <- tr$table$padj < 0.05
tr$table <- as.data.frame(tr$table) %>% bind_cols(feature = rownames(tr$table)) %>% mutate(log_padj = - log(.data$padj, base = 10))

FC = 0.6
colnames(tr$table) <- c("baseMean","logFC", "lfcSE","stat","PValue","FDR","threshold")
tr$table[tr$table$logFC > FC & tr$table$FDR < 0.05, "Expression"] <- "Up"
tr$table[tr$table$logFC < -FC & tr$table$FDR < 0.05, "Expression"] <- "Down"
tr$table[tr$table$FDR >= 0.05 , "Expression"] <- "Not Sig"
tr$table[tr$table$logFC > -FC & tr$table$logFC < FC, "Expression"] <- "Not Sig"

result$table <- tr$table[tr$table$threshold %in% T,]


y <- list("y"= y,"tr"= tr, "result"=result)


# Volcano plot
data <- y[['tr']]$table

FC=2
xscale = 5
yscale=4
seuil=0.6

if(length(unique(data$Expression))==2){color = c("grey","#F8766D")}else{color = c("#619CFF","grey","#F8766D")}
filter <- as.data.frame(data) %>% filter(-log10(data$FDR)>1.2 & abs(data$logFC)>FC) %>% arrange(FDR)

ggplot(data) +
  geom_point(aes(x = logFC, y = -log10(FDR), colour = Expression)) + 
  #ggrepel::geom_label_repel(data=filter[1:20,], aes(x = logFC, y = -log10(FDR))) +
  ggtitle("\nVolcano plot WT vs KO : FDR < 0.05 & logFC > 1.2\n") + xlab("log2 FC") + ylab("-log10 FDR") +
  scale_x_continuous(limits = c(-xscale,xscale)) + scale_y_continuous(limits = c(0,yscale)) + scale_color_manual(values = color) + 
  geom_hline(yintercept=1.30103, col="black", linetype="dashed") + geom_vline(xintercept=c(-seuil,seuil), col="black", linetype="dashed") +
  theme_classic() + theme(plot.title = element_text(size = rel(1.5), hjust = 0.5), axis.title = element_text(size = rel(1))) 

# PCA
rld <- rlog(y$y, blind = TRUE)
pca_data <- assay(rld)
pca <- prcomp(t(pca_data))
summary(pca)

library(ggplot2)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_df$condition <- colData(dds)$condition
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA des échantillons", x = "PC1", y = "PC2") +
  theme_minimal()


#3d
pca_df$PC3 <- pca$x[,3]
ggplot(pca_df, aes(x = PC1, y = PC2, z = PC3, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA 3D des échantillons", x = "PC1", y = "PC2", z = "PC3") +
  theme_minimal()
