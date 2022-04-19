
#Importing libraries
library(dada2)
library(ggplot2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(phangorn)
library(DESeq2)
#defining path
path = '.'
#importing merged raw reads
dataM <- sort(list.files(paste(path,'\\Raws',sep = ''), pattern=".fastq"))
#Extract sample names
sample.names <- sapply(strsplit(basename,(dataM), ".fastq"), `[`, 1)
sample.names
#Check quality profile of fastq files
plotQualityProfile(dataM[1:2])
#Creating filtered files
filt.data <- file.path(dataM, paste0(sample.names, "_filt.fastq"))
names(filt.data) <- sample.names
#filtrations
out <- filterAndTrim(dataM, filt.data, truncLen=c(420),
                      maxN=0, maxEE=c(2), truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE, truncLen value depends on quality of fastQ plot files
head(out)
#output filtration statestics
write.table(out,file=paste(path,'\\plots\\filtring_stat.tsv',sep=''),sep='\t')
#Learn error rates
err <- learnErrors(filt.data, multithread=TRUE)
#Assume the error rates are 3 times higher than what was estimated from the merged reads.
errM <- inflateErr(getErrors(err), 3)
#Plotting errors estimation
plotErrors(err, nominalQ=TRUE)
#Apply core sample inference algorithm on the trimmed reads
dada_merged <- dada(filt.data8, err=errM, multithread=TRUE)
#Construct amplicon sequence variant table
seqtab <- makeSequenceTable(dada_merged)
dim(seqtab8)
#Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#View dimensions
dim(seqtab.nochim)
#View percentage of chimeras removed
sum(seqtab.nochim)/sum(seqtab) * 100
#exporting table of ASVs without chimera
write.table(seqtab,paste(path,'\\results\\seqtab_nochim.tsv',sep=''),sep='\t')
#Track reads through the DADA2 pipeline
getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dada_merged, getN), rowSums(seqtab.nochim))
colnames(track.nbr.reads) <- c("input", "filtered", "denoised", "nonchim")
rownames(track.nbr.reads) <- sample.names
head(track.nbr.reads)
#Assigning taxonomy
taxa = assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fasta", multithread=TRUE)
#extracting sequences from ASvs
seqs = getSequences(seqtab.nochim)
names(seqs) = seqs
#Save taxonomy ASVs table in relevant directory
write.csv(taxa, file="~/results/ASVs_taxonomy.csv")
saveRDS(taxa, "~/results/ASVs_taxonomy.rds")
#Saving taxonomy ASVs as RDS object
tax <- readRDS("~/results/taxonomy.rds")
#Generating ASVs count table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
count.asv.tab <- t(seqtab.nochim)
head(count.asv.tab)
head(rownames(count.asv.tab))
#Saving ASVs count table
write.csv(count.asv.tab, file="~/results/ASVs_counts.csv")
saveRDS(count.asv.tab7, file="~/results/ASVs_counts.rds")
#Saving ASVs count table as RDS object
ASVs <- readRDS("~/results/ASVs_counts.rds")
#View dimensions of taxonomy table and ASV_counts
dim(tax)
dim(count.asv.tab) 
identical(rownames(count.asv.tab7),rownames(tax))  #TRUE

#Assign user-friendly ASV IDs to replace sequences
head(rownames(ASVs))
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
head(ASV.IDs)  #done!
#Now associate a given sequence to each pasted ASV using named vectors
names(seqs) <- ASV.IDs
head(seqs)
seq_lens <- nchar(seqs)#to recall the number of bps in each seq, we use *nchar* function
seq_lens
plot(density(seq_lens)) #check the sequence length graphically

#Merge ASV table and taxonomic table as a phyloseq object
phy <- phyloseq(otu_table(count.asv.tab7,taxa_are_rows = TRUE),tax_table(tax))
identical(taxa_names(phy),rownames(ASVs)) #TRUE
taxa_names(phy) <- names(seqs)
str(phy)
saveRDS(phy, "~/results/phy_rds")
tax_df<- data.frame(tax_table(phy))
head(tax_df)
dim(tax_df)
View(tax_df)
saveRDS(tax_df, "~/results/tax_ASVs.rds")


