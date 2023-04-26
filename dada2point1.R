#packages used
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.8")
library(dada2)
library(phyloseq)
library(csv)
library(tidyverse)


#-----MAKING ASV AND TAXA TABLES-----#

#path to fastq files
path <- "~/data/Oil_ML/Great_Lakes_Oil/MTUfastq/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

#fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S"), `[`, 1)
head(sample.names)

#visualize quality profiles
plotQualityProfile(fnFs[2:5])           #forward reads
plotQualityProfile(fnRs[2:5])           #reverse reads

#place filtered files in "filtered", a subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

#standard filtering parameters: 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 17, trimRight = 24,
                     maxN=0, maxEE=Inf, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


read.stats <- out %>%
  as_tibble() %>%
  summarise(
    max.in = max(reads.in),
    min.in = min(reads.in),
    mean.in = mean(reads.in),
    median.in = median(reads.in),
    sum.in = sum(reads.in),
    count.in = sum(ifelse(reads.in > 0, 1, 0)),
    under1000.in = sum(ifelse(reads.in < 1000, 1, 0)),
    max.out = max(reads.out),
    min.out = min(reads.out),
    mean.out = mean(reads.out),
    median.out = median(reads.out),
    sum.out = sum(reads.out),
    count.out = sum(ifelse(reads.out > 0, 1, 0)),
    under1000.out = sum(ifelse(reads.out < 1000, 1, 0))
  )
View(t(read.stats))

#save csvs
#write.csv(out, "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/dada2_output/out.csv")
#write.csv(read.stats, "/home/lgschaer/old/Plastic_Deg/DCPET_Zymo/samples2/dada2_output/read.stats.csv")

#learn about error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#inspecting dada-class object
dadaFs[[1]]
dadaRs[[1]]

#merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose=TRUE)

#inspect the merger data.frame from the first sample
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

View(track)


#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/data/sharedDatabases/Silva/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

view(taxa)


#save sequence table and taxa table
saveRDS(seqtab.nochim, "~/data/Oil_ML/Great_Lakes_Oil/MTUfastq/seqtab.rds")     #sequence table
saveRDS(taxa, "~/data/Oil_ML/Great_Lakes_Oil/MTUfastq/taxa.rds")                 #taxa table

seasonsSeqtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/Seasons/seqtab.rds") 
straitsSeqtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/EPA_fastq/seqtab.rds")
mtuSeqtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/MTUfastq/seqtab.rds")

seasonsTaxtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/Seasons/taxa.rds") 
straitsTaxtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/EPA_fastq/taxa.rds")
mtuTaxtab <- readRDS("~/data/Oil_ML/Great_Lakes_Oil/MTUfastq/taxa.rds")


meta <- read_csv("~/data/Oil_ML/Great_Lakes_Oil/SraRunInfo_EPA_MTU.csv", skip=1)
view(head(meta))

#meta <- data.frame(meta)

straitsSeqframe <- rownames_to_column(as.data.frame(straitsSeqtab))
mtuSeqframe <- rownames_to_column(as.data.frame(mtuSeqtab))

straitsSeqframe['Run'] <- as.data.frame(str_split(straitsSeqframe$rowname, "_", simplify = TRUE))$V1
mtuSeqframe['Run'] <- as.data.frame(str_split(mtuSeqframe$rowname, "_", simplify = TRUE))$V1


view(head(straitsSeqframe['Run']))
view(head(mtuSeqframe['Run']))

meta2 <- meta %>% select(c("Run", "SampleName"))


straitsFrame <- left_join(straitsSeqframe,meta2, by="Run")
mtuFrame <- left_join(mtuSeqframe,meta2, by="Run")

straitsFrame <- data.frame(straitsFrame, row.names="SampleName")
mtuFrame <- data.frame(mtuFrame, row.names="SampleName")
view(head(mtuFrame))

write.csv(straitsFrame, "~/data/Oil_ML/Great_Lakes_Oil/straitsFrame.csv")
write.csv(straitsTaxtab, "~/data/Oil_ML/Great_Lakes_Oil/straitsTaxtab.csv")

write.csv(seasonsSeqtab, "~/data/Oil_ML/Great_Lakes_Oil/seasonsFrame.csv")
write.csv(seasonsTaxtab, "~/data/Oil_ML/Great_Lakes_Oil/seasonsTaxtab.csv")

write.csv(mtuFrame, "~/data/Oil_ML/Great_Lakes_Oil/mtuFrame2.csv")
write.csv(mtuTaxtab, "~/data/Oil_ML/Great_Lakes_Oil/mtuTaxtab2.csv")