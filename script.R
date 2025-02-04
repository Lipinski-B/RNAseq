library(jsonlite)
library(dplyr)
library(ggplot2)

logs.path <- "I:/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/result/01_KALLISTO/"
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



















abundance <- "I:/ANALYSIS/BIT/PROJECTS/BL/Projet/SKILLS2/result/01_KALLISTO/Treated1_SE_abundance.tsv"
count <- read.table(file=abundance, sep = '\t', fill = T, header = T)
count <- count[count$tpm!=0,]
