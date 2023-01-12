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
library(dada2)
packageVersion("dada2")
#1.26.0
library(ShortRead)
packageVersion("ShortRead")
#1.56.1
library(Biostrings)
packageVersion("Biostrings")
#2.66.0

#make sure in correct directory
getwd()

#want to be in "/Users/samanthabeal/Documents/MSc/Bioinformatics"
setwd("Documents/MSc/Bioinformatics")

#directory should contain all sequence files (unzipped, primers removed)
path <- "dada2"
list.files(path)

#housekeeping
#generate matched lists of the forward and reverse read files, as well as parsing out the sample name
# Forward fastq file names have format: SAMPLENAME_L001_R1_001.fastq
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq", full.names = TRUE))

#identify primers
FWD <- "TAGCTCAGCTGGTAGAGCAC"
REV <- "TGCCATTTAGCAGACGCTTTT"

#check
#verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#pull out ambiguous bases -- not sure I need to but it's in the tutoral
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

#count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))

#cutadapt -- remove primers
#run 'whereis cutadapt' in terminal to get the path
cutadapt <- "/Library/Frameworks/Python.framework/Versions/3.10/bin/cutadapt"
system2(cutadapt, args = "--version")
#3.5

#If the above command succesfully executed, 
#R has found cutadapt and you are ready to continue following along.

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i])) # input files
}

#check
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]))

#i have 1 forward primer left somewhere
#from https://github.com/benjjneb/dada2/issues/675: a small number of primers left is ok to keep going w the subsequent analysis

#samples are now ready to be analyzed with dada2
# Forward and reverse fastq filenames have the format SAMPLENAME_L001_R1_001.fastq:
cutFs <- sort(list.files(path.cut, pattern="_L001_R1_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
#this code makes the sample names "190222-A1-1-42" (ie. removes everything after the first _)

#plot read quality profiles
plotQualityProfile(cutFs[1:2])
#quality drops off around 100bp (expected amplicon length = 90)
#will not be trimming to a fixed length

#Assigning the filenames for the output of the filtered reads to be stored as fastq.gz files.
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
#first tutorial had unzipped files, this tutorial is zipped.
#my files are unzipped so will need to run again to see if un/zipped makes a difference 

#Filter and trim
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 2, 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)
#file names are not shortened here but they are in the tutorial @_@

#learn error rate
errF <- learnErrors(filtFs, multithread = TRUE)
#101568798 total bases in 885772 reads from 23 samples will be used for learning the error rates.

#plot error rates
plotErrors(errF, nominalQ=TRUE)
#looks different then in the first tutorial (issues removing primers)
#A2A, C2C,& T2T points are not smooth

#Dereplicate identical reads
derepFs <- derepFastq(filtFs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#Sample Inference
#At this step, the core sample inference algorithm is applied to the dereplicated data.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

#Merge paired reads
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#skipping ^^

#Construct sequence table
#from https://github.com/benjjneb/dada2/issues/795: 
#if using single-end reads, make sequence table out of the dadaFs object, i.e. makeSequenceTable(dadaFs).
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
#51, 1674 (51 samples, 1674 ASVs)

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#59 identified and removed

#Sequence length distribution
table(nchar(getSequences(seqtab.nochim)))
#range from 51-151bp, majority are 151

#track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", 
                     "nonchim")
rownames(track) <- sample.names
head(track)








#######################
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
#101,419,010 total bases in 921,991 reads from 35 samples will be used for learning the error rates.

#visualize the estimated error rates
plotErrors(errF, nominalQ=TRUE)

#we are now ready to apply the core sample inference algorithm to the filtered and trimmed sequence data
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaFs[[1]]
#for sample 1 (190222-A1-1-42_F_filt.fastq.gz): 24 sequence variants were inferred from 13814 input unique sequences.
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
#SEQTAB.NOCHIM IS THE DATA I WANT! #READS/SAMPLE/ASV

sum(seqtab.nochim)/sum(seqtab)
#0.9951525
#~1% of sequence reads = chimeras

#track TOTAL reads throughout the pipeline
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
#phyloseq analysis----the dada2 tutiorial links to phyloseq for sequence analysis
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



##### assess seqtab.nochim data-----
#make a backup----
seqtab.nochim1 <- seqtab.nochim
seqtab.nochim1 <- as.data.frame(seqtab.nochim1)

#shorten names to be easier to look at -- this code keeps the first index (since cannot have duplicate row names)
rownames(seqtab.nochim1) <- stringr::str_extract(samples.out, "[^-]*-[^-]*-[^-]*")

#transpose so can see the read depth of each ASV/sample
t_seqtab.nochim <-  t(seqtab.nochim1)
t_seqtab.nochim <- as.data.frame(t_seqtab.nochim)


####working####
#need to translate the sequences into ASVs
#I want to blast the sequences, get the ASV info, and add it back into the file with the read counts
#separate out the sequences into a fasta file, can then Blast them
ASVs <- t_seqtab.nochim
ASVsfull <- ASVs 

library(tibble)
ASVs <- tibble::rownames_to_column(ASVs, "Seq")

ASVs <- ASVs[c(1)]
ASVs$Seq.name <- 1:nrow(ASVs) 

#switch order of columns
ASVs <- ASVs[,c(2,1)]

#seqs only
ASVseq <- ASVs[c(2)]
#make seqs into a list
ASVseqlist <- as.list(ASVseq)

install.packages("seqinr")
library("seqinr")

write.fasta(ASVseq, "SmaI corII", "dada2/dada2_SINE_gDNA.fasta", open = "w", nbchar = 110, as.string = FALSE)
#not quite

install.packages("phylotools")
library(phylotools)
dat2fasta(ASVs, outfile = "dada2/dada2_SINE_gDNA.fasta")


