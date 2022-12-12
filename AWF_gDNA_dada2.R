###########################
#dada2
###########################
#do this the first time
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install(version = '3.16')

#BiocManager::install("dada2", version = "3.16")
#packageVersion("dada2")
#1.26.0
###########################

#load libaries
library(BiocManager)
library(dada2)
library(Rcpp)

#make sure in correct directory
getwd()
#want to be in "/Users/samanthabeal/Documents/MSc/Bioinformatics"
setwd("Documents/MSc/Bioinformatics")

#directory should contain all sequence files
#get unzipped sequence files that have had the primers removed (taken from output of UNIX pipeline)
path <- "dada2"
list.files(path)

# Forward fastq filenames have format: SAMPLENAME_L001_R1.fastq and SAMPLENAME_L002_R1.fastq

fnFs <- sort(list.files(path, pattern="_L001_R1.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:2])

#quality drops off around 110
#make directory for filtered seqs (path = "Documents/MSc/Bioinformatics/dada2")
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

#trimming to 110 because that is where the sequence quality drops off as per the plot (approx.)
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
#24 sequence variants were inferred from 13814 input unique sequences.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16


#from here we would merge the filtered and denoised forward and reverse reads
###I have unpaired reads so will skip this step###

#mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

#make sequence table
#since did not merge, make table out of filtered R1 only
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
#51 samples, 487 ASVs

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#110,487 (487 ASVs all 110 bp long - because I trimmed to that length)


#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#Identified 21 bimeras out of 487 input sequences.
#51, 466 (51 samples, 466 ASVs)

sum(seqtab.nochim)/sum(seqtab)
#0.9951525
#~1% of sequence reads = chimeras

#track reads throughout the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

#assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "dada2/SINE_taxonomy.fasta", multithread=TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#check
track
taxa.print
#the track df tells me how many reads/sample, but not how many ASVs those reads belong to
#the taxa.print df tells me how many ASVs I have, but not the depth or which samples they came from

###########################
#phyloseq analysis----
###########################

library(phyloseq)
#packageVersion("phyloseq")
## [1] '1.42.0'
library(Biostrings)
#packageVersion("Biostrings")
## [1] '2.66.0'
library(ggplot2)
#packageVersion("ggplot2")
## [1] '3.4.0'
theme_set(theme_bw())

#make a simple dataframe out of the info in the sample name (ie. fish ID)
samples.out <- rownames(seqtab.nochim)
subject <- stringr::str_extract(samples.out, "[^-]*-[^-]*")
samdf <- data.frame(Fish.ID=subject)
rownames(samdf) <- samples.out


#construct a phyloseq object directly from the dada2 outputs.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#rename to be shorter
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

hist(sample_sums(ps), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=12)

#plot diversity (this is not really what I'm looking for but it's in the tutorial)
plot_richness(ps, x="Fish.ID", measures=c("Shannon", "Simpson"))

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Fish.ID")




