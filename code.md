## INTRODUCTION
Charger les librairies R utiliseés dans ce workshop
```
library(dada2)
library(phyloseq)
library(ggplot2)
```

Changer _YourPath_ par l'emplacement du dossier dada2_Analysis sur votre ordinateur
```
#MAC
output_directory = "YourPath/dada2_Analysis"
rawReadsFolder = sprintf("%s/rawReads",output_directory)
trainset_path = sprintf("%s/smaller_silva_train_set_128.fa.gz",output_directory)
designFile_path = sprintf("%s/Metadata/design.tsv",output_directory)

#WINDOWS
output_directory = "YourPath\\dada2_Analysis"
rawReadsFolder = sprintf("%s\\rawReads",output_directory)
trainset_path = sprintf("%s\\smaller_silva_train_set_128.fa.gz",output_directory)
designFile_path = sprintf("%s\\Metadata/design.tsv",output_directory)
```


## CHARGEMENT DES SEQUENCES
  - Détecter les fichiers contenant les séquences ADN
```
readFiles <- list.files(rawReadsFolder, pattern = "_R[12].fastq$", full.names = TRUE, recursive = TRUE)
length(readFiles)
#Record the sample names
sample.names <- unique(sapply(strsplit(basename(readFiles), "_R[12].fastq"), `[`, 1))
sample.names
```


  - Differencier les séquences _Forward_ des séquences _Reverse_
```
fnFs <- sort(list.files(rawReadsFolder, pattern = "_R1.fastq$", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(rawReadsFolder, pattern = "_R2.fastq$", full.names = TRUE, recursive = TRUE))
length(fnFs)
length(fnRs)
```

## PROFILS DE QUALITE
  - Examiner les profils de qualité de chacune des séquences _Forward_
```
dir.create(file.path(output_directory, subDir), showWarnings=FALSE)
savefolder = paste(output_directory, "/qualityProfilesRawReads", sep="")
#2. create the quality profile
cat("Forward reads\n")
for (i in fnFs) {
  j<-sapply(strsplit(i, sprintf("%s/",rawReadsFolder)), `[`, -1)
  destfile=sprintf("%s/qualityProfile_%s.pdf",savefolder,j)
  if (!file.exists(destfile)) {
    cat(".")
    pdf(file=destfile, width=10, height=7)
    print(plotQualityProfile(i))
    dev.off()
  }
}
```
  - Meme chose pour les séquences _Reverse_
```
for (i in fnRs) {
  j <- sapply(strsplit(i, sprintf("%s/",rawReadsFolder)), `[`, -1)
  destfile=sprintf("%s/qualityProfile_%s.pdf",savefolder,j)
  if (!file.exists(destfile)) {
    cat(".")
    pdf(file=destfile, width=10, height=7)
    print(plotQualityProfile(i))
    dev.off()
  }
}
```



## FILTRAGE DES SEQUENCES DE MAUVAISE QUALITE
  - Créer un dossier ou iront les sequences filtrees
  ```
filt_path <- file.path(output_directory, "Filtered_reads")
if(!file_test("-d", filt_path)) dir.create(filt_path)
```
  - Donner de nouveaux noms aux séquences filtrees
  ```
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
```
  - D'apres les observations des profils de qualite crees ci-dessus et en tenant compte de la longueur attendue des amplicons, donner des valeurs de coupe pour les sequences _Forward_ et _Reverse_
  ```
Fwd_trim=240
Rev_trim=220
```
  - Appliquer le filtre
  ```
#MAC/LINUX
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(Fwd_trim,Rev_trim),
                     maxN=0, compress=FALSE, multithread=TRUE) 
#WINDOWS
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(Fwd_trim,Rev_trim),
                     maxN=0, compress=FALSE, multithread=FALSE) 
```

## MACHINE LEARNING: auto-correction des sequences d'apres un modele parametrique
  - Creer un fichier pour le graphique genere ci-dessous
```
errorLearning_path <- file.path(output_directory, "Error_Rate_Learning")
if(!file_test("-d", errorLearning_path)) dir.create(errorLearning_path)
```
  - Lancer le modele parametrique
```
MAC:
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
WINDOWS:
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
```
  - Generer des graphiques de correction de sequence pour les sequence _Forward_ et _Reverse_
```
pdf(file=sprintf("%s/ErrorLearning_ForwardReads.pdf",errorLearning_path), width=10, height=7)
print(plotErrors(errF, nominalQ=TRUE))
dev.off()
pdf(file=sprintf("%s/ErrorLearning_ReverseReads.pdf",errorLearning_path), width=10, height=7)
print(plotErrors(errR, nominalQ=TRUE))
dev.off()
```



## DEREPLICATION
  - La dereplication va grouper les sequences identiques pour allerger la charge de calcul
```
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
```
  - Donner le nom des echantillons aux groupes de sequences derepliquees
```
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```



## INFERENCE DES ECHANTILLONS
  - Utiliser les resultat du modele parametrique precedent pour creer des sequences modifiees dans les echantillons  
```
#MAC/LINUX
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=FALSE)

#WINDOWS
dadaFs <- dada(derepFs, err=errF, multithread=FALSE, pool=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE, pool=FALSE)
```
  - Regarder combien de sequences ont ete creees dans chacun des echantillons
```
dadaFs[[1]]
dadaRs[[1]]
```

## ASSEMBLER LES SEQUENCES _FORWARD_ ET _REVERSE_ POUR OBTENIR LES AMPLICONS
  - Creer les amplicons
 ```
amplicons <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=20, verbose=TRUE)
 ```
  - Contruire un tableau avec  les sequences des amplicons
 ```
seqtab <- makeSequenceTable(amplicons)
 ```
  - Creer un histogram representant les differentes longueurs des amplicons
 ```
sequenceInference_path <- file.path(output_directory, "sequenceInference")
if(!file_test("-d", sequenceInference_path)) dir.create(sequenceInference_path)
pdf(file=sprintf("%s/Sequence_inference_histogram.pdf",sequenceInference_path), width=10, height=7)
hist(nchar(getSequences(seqtab)),  
     main=paste("Amplicon lengths overview:", dim(seqtab)[-1],"generated amplicons in total"),
     xlab="Amplicon length",
     border="black", col="gold",
     breaks=dim(seqtab)[-1])
dev.off()
 ```

# ENLEVER LES CHIMERES
 ```
#MAC/LINUX
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#WINDOWS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

```
  - Afficher le pourcentage de sequences detectees comme chimeres
```
cat("\nFraction of chimeras:", (1-sum(seqtab.nochim)/sum(seqtab))*100,"% of the total sequence reads.\n")
```


# DECOMPTE FINAL DES SEQUENCES
  - Creer un dossier sera creee le tableau de decompte des sequences
```
Summary_path <- file.path(output_directory, "Summary")
if(!file_test("-d", Summary_path)) dir.create(Summary_path)
```
  - Rassembler les statistiques de chacune des etapes
  ```
#Gather stats from the different steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
```
  - Creer le tableau de resultats
  ```
write.table(track ,file =sprintf("%s/SummaryTable.tsv",Summary_path), sep = "\t", quote = FALSE)
```

## ASSIGNATION DE LA TAXONOMIE
```
taxa <- assignTaxonomy(seqtab.nochim, trainset_path)
```

## TABLEAU DE COMPTE DES AMPLICONS ET TABLEAU D'IDENTIFICATION TAXONOMIQUE
- Creer un dossier seront crees ces fichier
```
dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
if(!file_test("-d", dada2outputFiles_path)) dir.create(dada2outputFiles_path)
```
  - Tableau de compte des amplicons
  ```
write.table(t(seqtab.nochim), file =sprintf("%s/otu_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)
```
- Tableau d'identification taxonomique
```
write.table(taxa ,file =sprintf("%s/tax_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)
```



## EXPLORATION DES DONNEES AVEC PHYLOSEQ
  - Creer un dossier pour les graphiques
  ```
graphs_path <- file.path(output_directory, "graphics")
if(!file_test("-d", graphs_path)) dir.create(graphs_path)
```
  - Integrer les donnees de Dada2 et du fichier design pour creer un object Phyloseq
1. Design de l'experience
```
presampleTable = read.table(file =designFile_path, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
presampleTable<- as.data.frame(presampleTable)
sampleTable = sample_data(presampleTable)
```
2. Creer un object otu_table pour phyloseq
```
OtuTable = otu_table(seqtab.nochim, taxa_are_rows=FALSE)
```
3. Creer un object tax_table pour phyloseq
```
TaxTable = tax_table(taxa)
```
4. Creer l'objet phyloseq
```
ps <- phyloseq(OtuTable, TaxTable, sampleTable)
ps
```

## Diversite alpha (abondance des especes)
  - Rarefaction: pour comparer plusieurs echantillons ensemble
 ```
psRar = rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),verbose = TRUE)
```

- Extraire les 50 especes les plus abondantes
```
top50 <- names(sort(taxa_sums(psRar), decreasing=TRUE))[1:50]
psRar.top50 <- transform_sample_counts(psRar, function(OTU) OTU/sum(OTU))
psRar.top50 <- prune_taxa(top50, psRar.top50)
listTaxLevels=c("Phylum", "Class", "Order", "Family", "Genus", "Species")
destfile=sprintf("%s/barCharts_top50.pdf",graphs_path)
pdf(file=destfile, width=12, height=9)
for(l in listTaxLevels) {
  p = plot_bar(psRar.top50, x=sprintf("%s",l), fill=sprintf("%s",l), facet_grid=~Conditions, title=sprintf("%s level of the 20 most abundant ASVs", l))
  q= p + geom_bar(aes_string(color=sprintf("%s",l), fill=sprintf("%s",l)), stat="identity", position="stack")
  print(q)
}
dev.off()
```




## Diversite Alpha (Richesse spécifique)
```
destfile=sprintf("%s/richness.pdf",graphs_path)
pdf(file=destfile, width=12, height=12)
plot_richness(ps, x = "Conditions", color = "Conditions", measures=c("Observed", "Chao1", "Shannon", "Simpson"), nrow=1) + 
  geom_boxplot() + geom_point(size = 0.5, alpha = 0.5) +
  ggtitle("Richness plot") +
  theme_classic(base_size = 16)
dev.off()
```


## Diversite Beta (ordination): 3 mesures differentes
```
#Fix for one condition experiments
sample_data(psRar)[ , 2] <- sample_data(psRar)[ ,1]
destfile=sprintf("%s/ordination.pdf",graphs_path)
pdf(file=destfile, width=10, height=7)
print("PCoA")
psRar.pcoa <- ordinate(psRar, method="PCoA", distance="bray")
p=plot_ordination(psRar, psRar.pcoa, type="samples", color="Conditions")
df=p$data
ggplot(df, aes(df[,1], df[,2], color = Conditions)) + ggtitle("PCoA plot") + xlab(colnames(df[1])) + ylab(colnames(df[2])) +
  geom_point(size = 3, alpha = 0.8) +geom_text(aes(label=row.names(df)),hjust=0, vjust=0, size = 3, color = 'grey') +  theme_classic(base_size = 16)
print("CCA")
psRar.cca <- ordinate(psRar, method="CCA", distance="bray")
p=plot_ordination(psRar, psRar.cca, type="Sample", color="Conditions")+
  ggtitle("CCA")
df=p$data
ggplot(df, aes(df[,1], df[,2], color = Conditions)) + 
  ggtitle("CCA plot") + xlab(colnames(df[1])) + ylab(colnames(df[2])) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label=row.names(df)),hjust=0, vjust=0, size = 3, color = 'grey') +
  theme_classic(base_size = 16)
print("MDS")
temp = psRar
otu_table(temp)[otu_table(temp) < 0.0] <- 0.0
psRar.mds <- ordinate(temp, method="MDS", distance="bray")
p=plot_ordination(temp, psRar.mds, type="Sample", color="Conditions")+
  ggtitle("MDS")
df=p$data
ggplot(df, aes(df[,1], df[,2], color = Conditions)) + 
  ggtitle("MDS plot") + xlab(colnames(df[1])) + ylab(colnames(df[2])) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label=row.names(df)),hjust=0, vjust=0, size = 3, color = 'grey') +
  theme_classic(base_size = 16)
dev.off()
```
