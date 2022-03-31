
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
#importing raw files
fnFs <- sort(list.files(paste(path,'\\Raws',sep = ''), pattern="_R1_001.fastq"))
fnRs <- sort(list.files(paste(path,'\\Raws',sep = ''), pattern="_R2_001.fastq"))
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
#Creating filtred files
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))
#filtrations
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
#output filtration statestics
write.table(out,file=paste(path,'\\plots\\filtring_stat.tsv',sep=''),sep='\t')
#estimating errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#ploting errors estimation
plotErrors(errF, nominalQ=TRUE)
#dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
#merging
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#exporting statestics of merging 
write.table(mergers[[1]],file = paste(path,'\\plots\\mergine_stat.tsv',sep=''),sep='\t')
#creating ASVs table
seqtab <- makeSequenceTable(mergers)
#Exporting ASVs table
write.table(seqtab,paste(path,'\\results\\seqtab.tsv',sep=''),sep='\t')
#removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#exporting table of ASVs without chimera
write.table(seqtab,paste(path,'\\results\\seqtab_nochim.tsv',sep=''),sep='\t')
#Assigning taxonomy
taxa = assignTaxonomy(seqtab.nochim, "silva_nr_v132_train_set.fasta", multithread=TRUE)
#extracting sequences from ASvs
seqs = getSequences(seqtab.nochim)
names(seqs) = seqs
#Alignment
alignment = AlignSeqs(DNAStringSet(seqs), anchor=NA)
phang.align = phyDat(as(alignment, "matrix"), type="DNA") 
#constructing ML model
dm = dist.ml(phang.align)
treeNJ = NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR = update(fit, k=4, inv=0.2)
fitGTR = optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
            rearrangement = "stochastic", control = pml.control(trace = 0))
#Constructing Phyloseq object
ps = phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa),phy_tree(fitGTR$tree))
#ploting diversity 
plot_richness(ps, x="Status", measures=c("Shannon", "Simpson"), color="Status") + theme_bw()
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
plot_ordination(ps, ord.nmds.bray, color="Status", title="Bray NMDS")
#Constructing DeSeq object
ds = phyloseq_to_deseq2(ps, ~ Status)
ds = DESeq(ds)
alpha = 0.01
res = results(ds, contrast=c("Status", "Case", "Control"), alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res_sig = res[(res$padj < alpha), ]
res_sig = cbind(as(res_sig, "data.frame"), as(tax_table(ps_nona)[rownames(res_sig), ], "matrix"))
ggplot(res_sig, aes(x=Genus, y=log2FoldChange, color=Family)) + geom_jitter(size=3, width = 0.2) + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
