###########################
#dada2
###########################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = '3.16')

BiocManager::install("dada2", version = "3.16")

library(dada2)
packageVersion("dada2")
#1.26.0


getwd()
setwd("Documents/MSc/Bioinformatics")

#get unzipped sequence files that have had the primers removed (taken from output of UNIX pipeline)
path <- "AWF_gDNA/dada2"
list.files(path)

# Forward fastq filenames have format: SAMPLENAME_L001_R1.fastq and SAMPLENAME_L002_R1.fastq

fnFs <- sort(list.files(path, pattern="_L001_R1.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])

#quality drops off around 110
#make directory for filtered seqs
filtFs <- file.path(path, "AWF_gDNA/dada2/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, truncLen=c(110),
                     maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#learn error rate
errF <- learnErrors(filtFs, multithread=TRUE)
#101419010 total bases in 921991 reads from 35 samples will be used for learning the error rates.

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

#we are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaFs[[1]]

#from here we would merge the filtered and denoised forward and reverse reads
#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

#I have unpaired reads so will skip this step
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
#51, 487 (487 ASVs)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#110,487 (487 ASVs all 110 bp long - because I trimmed to that length)


#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#51, 466

sum(seqtab.nochim)/sum(seqtab)
#0.9951525

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "AWF_gDNA/dada2/SINE_taxonomy.fasta", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#make into dataframe and combine into one file
track <- as.data.frame(track)
taxa.print <- as.data.frame(taxa.print)

#the track df tells me how many reads/sample, but not how many ASVs those reads belong to
#the taxa.print df tells me how many ASVs I have, but not the depth or which samples they came from

###########################
#phyloseq analysis----
###########################

library(phyloseq); packageVersion("phyloseq")
## [1] '1.42.0'
library(Biostrings); packageVersion("Biostrings")
## [1] '2.66.0'
library(ggplot2); packageVersion("ggplot2")
## [1] '3.4.0'
theme_set(theme_bw())

#construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))

#rename to be shorter
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps







