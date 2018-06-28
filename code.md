Charger les librairies R utilisees dans ce workshop
```
library(dada2)
library(phyloseq)
library(ggplot2)
```

Changer _YourPath_ par le chemin du dossier dada2_Analysis sur votre ordinateur
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


#Detect how many reads are in the rawReadsFolder
readFiles <- list.files(rawReadsFolder, pattern = "_R[12].fastq$", full.names = TRUE, recursive = TRUE)
length(readFiles)
#Record the sample names
sample.names <- unique(sapply(strsplit(basename(readFiles), "_R[12].fastq"), `[`, 1))
sample.names

# Forward and reverse fastq filenames have the following format: SampleName.pair[12].fastq
fnFs <- sort(list.files(rawReadsFolder, pattern = "_R1.fastq$", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(rawReadsFolder, pattern = "_R2.fastq$", full.names = TRUE, recursive = TRUE))
length(fnFs)
length(fnRs)


#Examine quality profiles graphs of forward and reverse reads
#1. Create a qulaity profile folder
subDir = sprintf("qualityProfilesRawReads")
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
cat("\nReverse reads\n")
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
cat("\n-----\ndone! Check pdf files in: ", savefolder,"\n", sep = "")




#FILTERING STEP (remove low quality reads)
#1. create a folder for the filtered reads
filt_path <- file.path(output_directory, "Filtered_reads")
if(!file_test("-d", filt_path)) dir.create(filt_path)
#2. create names for the filtered reads
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
#Apply the filter
Fwd_trim=240 #from quality profiles observation + expected amplicon length
Rev_trim=220 #from quality profiles observation + expected amplicon length

#MAC/LINUX
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(Fwd_trim,Rev_trim),
                     maxN=0, compress=FALSE, multithread=TRUE) 
#WINDOWS
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(Fwd_trim,Rev_trim),
                     maxN=0, compress=FALSE, multithread=FALSE) 

#ERROR RATE LEARNING: PARAMETRIC ERROR MODEL AUTO-CORRECTION (MACHINE LEARNING)
#Error rates are learned by alternating between sample inference and error rate estimation until convergence
#1. Create a folder
errorLearning_path <- file.path(output_directory, "Error_Rate_Learning")
if(!file_test("-d", errorLearning_path)) dir.create(errorLearning_path)
#2. Apply the statistical model
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
#3. Create error rates graphs
pdf(file=sprintf("%s/ErrorLearning_ForwardReads.pdf",errorLearning_path), width=10, height=7)
print(plotErrors(errF, nominalQ=TRUE))
dev.off()
pdf(file=sprintf("%s/ErrorLearning_ReverseReads.pdf",errorLearning_path), width=10, height=7)
print(plotErrors(errR, nominalQ=TRUE))
dev.off()



cat("\n########################### DEREPLICATION #################################\n")

start_time <- Sys.time()
cat("\nStarting dereplication  process")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names




# DEREPLICATION 
#Combine all identical sequencing reads into into “unique sequences”
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the dereplicated objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


# SAMPLE INFERENCE
#Apply the core sequence-variant inference algorithm to the dereplicated data
#MAC/LINUX
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=FALSE)

#WINDOWS
dadaFs <- dada(derepFs, err=errF, multithread=FALSE, pool=FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=FALSE, pool=FALSE)
dadaFs[[1]]
dadaRs[[1]]


# OVERLAP PAIRED-END READS TO GETR THE AMPLICON
#Create amplicons
amplicons <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap=20, verbose=TRUE)
#2. Constructing sequence table")
seqtab <- makeSequenceTable(amplicons)
#3 Create amplicon length histogram
sequenceInference_path <- file.path(output_directory, "sequenceInference")
if(!file_test("-d", sequenceInference_path)) dir.create(sequenceInference_path)
pdf(file=sprintf("%s/Sequence_inference_histogram.pdf",sequenceInference_path), width=10, height=7)
hist(nchar(getSequences(seqtab)),  
     main=paste("Amplicon lengths overview:", dim(seqtab)[-1],"generated amplicons in total"),
     xlab="Amplicon length",
     border="black", col="gold",
     breaks=dim(seqtab)[-1])
dev.off()

# REMOVE TECHNICAL CHIMERAS
#MAC/LINUX
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#WINDOWS
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)


cat("\nFraction of chimeras:", (1-sum(seqtab.nochim)/sum(seqtab))*100,"% of the total sequence reads.\n")


write.table(seqtab.nochim, file =sprintf("%s/seqtab.nochim.txt",sequenceInference_path), sep = "\t", quote = FALSE)


# SUMMARY TABLE #################################\n")
#1. Create folder for the table
Summary_path <- file.path(output_directory, "Summary")
if(!file_test("-d", Summary_path)) dir.create(Summary_path)
#Gather stats from the different steps
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(amplicons, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#output the table
write.table(track ,file =sprintf("%s/SummaryTable.tsv",Summary_path), sep = "\t", quote = FALSE)


# TAXONOMY ASSIGNMENT
taxa <- assignTaxonomy(seqtab.nochim, trainset_path)


# CREATE MAIN OUTPUT: AMPLICON COUNT TABLE AND TAXONOMY TABLE
#1. Create a folder for these files
dada2outputFiles_path <- file.path(output_directory, "outputFilesDada2")
if(!file_test("-d", dada2outputFiles_path)) dir.create(dada2outputFiles_path)
#2 output amplicon count table
write.table(t(seqtab.nochim), file =sprintf("%s/otu_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)
#3 output taxonomy table
write.table(taxa ,file =sprintf("%s/tax_table.txt",dada2outputFiles_path), sep = "\t", quote = FALSE)




############################# DATA EXPLORATION   ##########################################
library(phyloseq)


# SOME GRAPHS
#1. create a folder for graphics
graphs_path <- file.path(output_directory, "graphics")
if(!file_test("-d", graphs_path)) dir.create(graphs_path)
#2. Load dada2 output into phyloseq by merging the 2 output and the design file into a same object
#  a. Load the design file as a table
presampleTable = read.table(file =designFile_path, header=T, sep = "\t", row.names=1, com='', check.names=FALSE)
#  b. convert it as a dataframe
presampleTable<- as.data.frame(presampleTable)
#  c. Create a sample_data-class object for phyloseq
sampleTable = sample_data(presampleTable)
#  d. Create a  otu_table-class object for phyloseq
OtuTable = otu_table(seqtab.nochim, taxa_are_rows=FALSE)
#  e. Create a taxonomyTable-class object for phyloseq
TaxTable = tax_table(taxa)
#3. Loading a phyloseq object
ps <- phyloseq(OtuTable, TaxTable, sampleTable)
ps


# Alpha diversity plots (species abundance)
#Rarefaction: comparing several samples together
psRar = rarefy_even_depth(ps, sample.size = min(sample_sums(ps)),verbose = TRUE)


#1. take the top 20 most abundant species
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





# Alpha diversity: Richness plots 
destfile=sprintf("%s/richness.pdf",graphs_path)
pdf(file=destfile, width=12, height=12)
plot_richness(ps, x = "Conditions", color = "Conditions", measures=c("Observed", "Chao1", "Shannon", "Simpson"), nrow=1) + 
  geom_boxplot() + geom_point(size = 0.5, alpha = 0.5) +
  ggtitle("Richness plot") +
  theme_classic(base_size = 16)
dev.off()



# Beta-diversity: Ordination plots
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

