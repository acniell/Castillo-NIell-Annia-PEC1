library(readr)
human_cachexia <- read_csv("2024-Cachexia/human_cachexia.csv")
View(human_cachexia)
dim(human_cachexia)
head(human_cachexia)

library(GEOquery)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
BiocManager::install("mixOmics")
library(SummarizedExperiment)
library(Biobase)
x<-as.matrix(human_cachexia[, -c(1, 2)]) 
columnesMeta <- human_cachexia[, 1:2] 
colnames(columnesMeta) <- c("ID pacient", "Muscle loss")
rownames(x) <- columnesMeta$`ID pacient` 
SumExp <- SummarizedExperiment(
  assays = list(counts = t(x)),  
  colData = DataFrame(columnesMeta), 
  rowData = DataFrame(Variable = colnames(x))  
)
SumExp
dataset<-assay(SumExp) 

colData(SumExp)
rowData(SumExp)
dim(SumExp)
head(colData(SumExp))
any(is.na(dataset))
summary(dataset)
sd(dataset)
summary(t(dataset))
sd(t(dataset))
par(mfrow=c(2, 4))  
for (i in 1:63) {  
  boxplot(x[, i] ~ columnesMeta$`Muscle loss`,
          main = colnames(x)[i],  
          xlab = "Pèrdua muscular",
          ylab = "Expressió",
          col = c("cadetblue3", "coral"),
          border ="darkgreen")
}
par(mfrow=c(1, 1)) 

plot(density(dataset[, 1]), main = "Metabòlits", xlab = "Quantitat de metabòlit", ylab = "Densitat de probabilitat", col = "blue", lwd = 2)

for (i in 2:ncol(dataset)) {
  lines(density(dataset[, i]), col = i, lwd = 2)  }

datasetLog2<- log2(assays(SumExp)$counts + 1)

plot(density(datasetLog2[, 1]), main = "Metabòlits (log2)", xlab = "Expressió (log2)", ylab = "Densitat", col = "blue", lwd = 2)

for (i in 2:ncol(datasetLog2)) {
  lines(density(datasetLog2[, i]), col = i, lwd = 2) 
}

log2_dataset_normalitzat <- scale(datasetLog2, center = TRUE, scale = TRUE)

PCA_calcul <- prcomp(t(log2_dataset_normalitzat), center = TRUE, scale. = TRUE)
summary(PCA_calcul)
PCA_calcul$x
library(mixOmics)
PCA_tune<-tune.pca(t(log2_dataset_normalitzat), ncomp=10, scale=FALSE)
plot(PCA_tune)
coord <- PCA_calcul$x
PCA_resultat <- as.data.frame(coord)
PCA_resultat$Muscle_loss <- colData(SumExp)$Muscle.loss
suppressMessages({
  ggplot(PCA_resultat, aes(x = PC1, y = PC2, color = Muscle_loss)) +
    geom_point(shape = 17) +
    theme_minimal() +
    labs(title = "Anàlisi de Components Principals", x = "PC1", y = "PC2")
})"Anàlisi de Components Principals", x = "PC1", y = "PC2")


t_test_PC1 <- t.test(PC1 ~ Muscle_loss, data = PCA_resultat)
print(t_test_PC1)
t_test_PC2 <- t.test(PC2 ~ Muscle_loss, data = PCA_resultat)
print(t_test_PC2)
loads <- PCA_calcul$rotation
plot(loads[, 1], loads[, 2],
     xlab = "Càrregues PC1",
     ylab = "Càrregues PC2",
     main = "Càrregues PC1 i PC2",
     pch = 19, col = "purple")
text(loads[, 1], loads[, 2], labels = rownames(loads), pos = 4, cex = 0.5)
loads<-PCA_calcul$rotation
PC2_loads <- sort(abs(loads[, "PC2"]), decreasing = TRUE)
metabolits_PC2 <- names(PC2_loads[1:10]) 
data.frame(Metabolit = metabolits_PC2, Carrega_PC2 = loads[metabolits_PC2, "PC2"])
PC2_dataframe <- data.frame(Metabolit = metabolits_PC2, Carrega = loads[metabolits_PC2, "PC2"])
ggplot(PC2_dataframe, aes(x = reorder(Metabolit, Carrega), y = Carrega)) +
  geom_bar(stat = "identity", fill = "orange") +
  coord_flip() +
  labs(title = "Càrregues dels metabòlits a PC2",
       x = "Metabolit",
       y = "Càrrega a PC2") +
  theme_minimal()
metabolit_data <- data.frame(Muscle_loss = columnesMeta$`Muscle loss`)
for (metabolit in metabolits_PC2) {
  metabolit_data[[metabolit]] <- dataset[metabolit, ]
}
results <- lapply(metabolits_PC2, function(metabolit) {
  t.test(metabolit_data[[metabolit]] ~ metabolit_data$Muscle_loss, var.equal = FALSE)
})
for (i in 1:length(metabolits_PC2)) {
  cat(metabolits_PC2[i], ":\n")
  print(results[[i]])
  cat("\n")
}

#### Gràfic 1: boxplot comparatiu de cada metabòlit entre caquèctic i control {#gràfic-1-boxplot-comparatiu-de-cada-metabòlit-entre-caquèctic-i-control}
#### Gràfic 2: Gràfic de densitat i distribució dels metabòlits {#gràfic-2-gràfic-de-densitat-i-distribució-dels-metabòlits}
#### Gràfic 3: Gràfic de densitat i distribució dels metabòlits amb la transformació logarítmica {#gràfic-3-gràfic-de-densitat-i-distribució-dels-metabòlits-amb-la-transformació-logarítmica}
#### Gràfic 4: Variància explicada per cada PC {#gràfic-4-variància-explicada-per-cada-PC}
#### Gràfic 5: Representació gràfica de PC1 i PC2 {#gràfic-5-representació-gràfica-de-PC1-i-PC2}
#### Gràfic 6: Càrregues de PC1 i PC2 {#gràfic-6-càrregues-de-PC1-i-PC2}
#### Gràfic 7: Càrregues dels metabòlits a PC2 {#gràfic-7-càrregues-dels-metabòlits-a-PC2}


