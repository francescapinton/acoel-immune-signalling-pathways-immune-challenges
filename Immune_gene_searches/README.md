# acoel-immune-signalling-pathways-immune-challenges

Used on Ubuntu 24.04.1 LTS

Python:  
Python 3.8.10  

HMMER:  
HMMER3 - Version: 3.3 (Nov 2019)  

Mafft:  
v7.511 (2022/Dec/15)

iqtree2:  
IQ-TREE multicore version 2.2.0 COVID-edition for Linux 64-bit built Jun  1 2022 

AMAS:
https://github.com/marekborowiec/AMAS
Borowiec, M.L. 2016. AMAS: a fast tool for alignment manipulation and computing of summary statistics. PeerJ 4:e1660.

## 01 From gene sequences to alignments of domains from CDD  
Input file: a Conserved Domain Batch search output downloaded file for sequences of the same gene of interest belonging to different species.  
It extracts domain names and accession numbers present in the file, then it asks the user to manually remove the unwanted ones.  
It then downloads alignments for each domain from Conserved Domain Database (CDD)  

## 02 Build hmm profile, search them against proteome and extract domains or sequences
Input: directory with alignments to search, proteomes  
It builds hmm profiles from alignments for each domain (hmmbuild by hmmer)  
It searches hmm profiles against the proteome (translated transcriptome or genome) (hmmsearch by hmmer)  
It extracts domain sequences from the proteome, based on envelope coordinates  
Output: Fasta files named after the domain, containing domain sequences with defline as in the proteome  

## 02bis Search hmm profiles against fasta files with multiple gene sequences and extract domains or sequences
Input: hmm profiles folder, fasta files with multiple sequences  
To search already built hmm profiles (from script 02 for example) against fasta files with multiple gene sequences  
It extracts domains and append them to files of interest  

As 02 but but without hmmbuild  
Useful to add domains from extra databases or if hmm profiles are already available

## 03 alignment, model selection and gene tree partitioned by domain
Input: files named after a domain, containing domain sequences of genes from various species, named after the gene (or sequence name in the transcriptome)  
It aligns sequences with MAFFT  
It concatenates the alignments and creates a partition file with AMAS  
It runs iqtree2 to perform best partition and best model search and build ML tree

## Extract_domains_or_seq_only_from_domtable  
Runs extract_domains.py to extract all domains or full sequences from directories ending in _hmm_search in the current directory.
