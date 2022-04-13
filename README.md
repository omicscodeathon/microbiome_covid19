# Covid19_Microbiome
Using a data mining approach to analyze COVID-19 microbime sequences.

## Background :dna:: 
In the year 2020, SARS-CoV-2 (coronavirus 2019, COVID-19) has created an exceptional circumstance. Despite the fact that numerous other viruses had surfaced in the past, none of them had caused a global emergency like this one. According to official sources, this virus started in China and has been rapidly spreading throughout the world. On March 11, the World Health Organization (WHO) proclaimed the pandemic, and a huge number of countries have had to put in place measures to combat infection transmission.

Human microbiotas are symbiotic populations of microorganisms that live with humans. They play a key function in the immunological response of the host to respiratory viral infection. However,there is inadequate evidence to link the human microbiome to coronavirus disease (COVID-19). Also, the disparities in COVID-19 microbiome composition between people of different demographic and clinical status – age, race, sex, gender, varying symptoms and countries remain uncertain at the moment

Antibiotic treatment is the primary therapeutic method utilized to treat COVID-19 infection. Given that such method quickly produces antibiotic-resistant strains of opportunistic microorganisms, improved antibiotic therapy are essential to effectively control long-term symptoms and future pandemics, notably in patients infected with SARS-CoV-2.

Playing a key role in precision medicine, the gut microbiome is one important component to study among different populations in order to develop effective adjuncent treatements for instance targeted probiotics or FMT that can be of great use for the COVID-19 infection or for long-term COVID one helping to ceize symptoms such as diarrhea. 


## Objective:
In this study we implemented a data mining and meta-analysis approach to cluster and analyze COVID-19 microbiome sequences from cross-study datasets. The goal was to compare the community structure across these datasets and to cluster the resulting microorganisms that were mostly abundant in or common between each country’s samples (USA, Italy, Germany, France, China, India, Japan, North America) and to see if some these microorganisms are linked to diseases. 

## Dependencies 🖥️:
*Software:*

[DADA2](https://benjjneb.github.io/dada2/): **Fast, accurate, single-nucleotide resolution for amplicon data**
The DADA2 algorithm for the inference of exact amplicon sequence variants (ASVs) from amplicon data is implemented in the dada2 R package available in Bioconductor. The core algorithm replaces the traditional OTU picking step in 16S/18S/ITS marker-gene surveys with the inference of the exact sequences present in the sample after errors are removed.

### Installation:
Binaries for the current release version of DADA2 (1.14) are available from Bioconductor. Note that you must have R 3.6.1 or newer, and [Bioconductor version 3.10](https://www.bioconductor.org/install/), to install the current release from Bioconductor.

`if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")`

## DADA2 workflow:
![worklow](https://github.com/omicscodeathon/microbiome_covid19/blob/main/figures/DADA2-Workflow.png) 

## Output:
**1- A feature table of amplicon sequence variants (an ASV table):**
A matrix with rows corresponding to samples and columns to ASVs, in which the value of each entry is the number of times that ASV was observed in that sample. 

**2- Taxonomy table:** 
Assign taxonomy to phylogenetically informative marker-gene data, such as the 16S or 18S rRNA gene and the ITS region in fungi. make taxonomic assignments by comparing to a set of taxonomically assigned sequences in a provided reference fasta file. Appropriately formatted reference fasta files for several popular reference databases are available, including Silva, GreenGenes, RDP and UNITE.

## Team members:
- [Sofia Sehli](https://github.com/SofSei). Mohammed VI University of Health Sciences of Casablanca, Morocco.
- [Nihal Habib](https://github.com/NihalHB). Mohammed VI University of Health Sciences of Casablanca, Morocco.
- [Adijat Jimoh](https://github.com/adijatj). Division of Immunology, Faculty of Health Sciences, Institute of Infectious Diseases and Molecular Medicine, University of Cape Town, South Africa.
 Department of Genetics, Genomics and Bioinformatics, National Biotechnology Development Agency, Abuja, Nigeria
