#!/usr/bin/env Rscript


## Takes BLAST (format 6) results from BLASTing against IMG/VR and
## resolves, where possible, the consensus taxonomy string for
## each query

## $1: Format 6 BLAST results. Tab-delimited.
## $2: GFF file from PRODIGAL for BLASTed sequences
## $3: Reference taxonomy file. Tab-delimited. Two columns : first ID, second corresponding taxonomy
## $4: Output file name

args = commandArgs(trailingOnly=TRUE)

library(tidyverse)

# blast_results.df <- read.table(args[1], sep="\t", header=F)
# imgvr_taxonomy.df <- read.delim("args[2], sep="\t", header=F)

# Load BLAST data
setwd("/Users/julianzaugg/Desktop/ACE/major_projects/CF_phage/minion/imgvr")
blast_results.df <- read.table("viral_proteins_imgvr_diamond_blast.tsv", sep = "\t", header = T)

# Load GFF file from prodigal
gff.df <- read.delim("all_samples_viral_sequences.gff", header=F, comment.char="#")
names(gff.df) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# Add protein ID to seq_id
gff.df$complete_seq_id <- with(gff.df, paste0(seqid,"_", gsub("ID=[0-9]{1,10}_(.*?);.*", "\\1", gff.df$attributes)))

# Load taxonomy reference file.
imgvr_taxonomy.df <- read.delim("IMGVR_taxonomy.tsv", header = T, sep = "\t")
rownames(imgvr_taxonomy.df) <- imgvr_taxonomy.df$UViG
# imgvr_taxonomy.df[grepl("Gammaproteobacteria_gi_444434353", imgvr_taxonomy.df$UViG),]

# Count the number of genes per sequence
Sequence_gene_count.df <-
  gff.df %>%
  dplyr::group_by(seqid) %>%
  dplyr::summarise(Number_of_Genes = n())

# Subject - lineage
# Clean Subject_ID to be able to search against taxonomy reference
# Everything after first '|' assumed to be specific reference protein ID
blast_results.df$Subject_ID_cleaned <- gsub("\\|.*","", blast_results.df$Subject_ID)

subject_ids.v <- unique(blast_results.df$Subject_ID_cleaned)
subjects_imgvr_taxonomy.df <- imgvr_taxonomy.df[subject_ids,]
subjects_imgvr_taxonomy.df <- subjects_imgvr_taxonomy.df[!is.na(subjects_imgvr_taxonomy.df$UViG),]

missing_subject_ids.v <- subject_ids.v[!subject_ids.v %in% rownames(imgvr_taxonomy.df)]

# Even after basic cleaning, some entries won't match exactly. Need to perform more robust searching
# Search for any match to the subject ID
# blast_results.df$Subject_ID_cleaned[!blast_results.df$Subject_ID_cleaned %in% imgvr_taxonomy.df$UViG]
# imgvr_taxonomy.df <- rbind(imgvr_taxonomy.df, c("Gammaproteobacteria_gi_444434353.6151-50330x", "xxx"))
# imgvr_taxonomy.df <- rbind(imgvr_taxonomy.df, c("Gammaproteobacteria_gi_444434353.6151-50330xx", "xxx"))
# imgvr_taxonomy.df[grepl("NC_003278[\\|\\.]", imgvr_taxonomy.df$UViG),]

missing_subject_matches.l <- list()
# output[[tax_string_level]] <- list("counts" = counts.m, "abundances" = abundances.m)
# print(list(grep(paste0("Gammaproteobacteria_gi_444434353", "[\\|\\.]"), imgvr_taxonomy.df$UViG, value =T)))
for (sid in missing_subject_ids.v){
  # print(grep(paste0(sid, "[\\|\\.]"), imgvr_taxonomy.df$UViG, value =T))
  missing_subject_matches.l[sid] <- list(grep(paste0(sid, "[\\|\\.]"), imgvr_taxonomy.df$UViG, value =T))
}



lapply(blast_results.df$Subject_ID_cleaned, function(x) grepl(blast_results.df$Subject_ID_cleaned))

# Sequence ID - Gene ID - lineage

temp[!temp %in% imgvr_taxonomy.df$UViG]
gsub("\\|.*","", blast_results.df$Subject_ID)  %in% imgvr_taxonomy.df$UViG

# blast_results.df$Percent_query_aligned <- round(abs(with(blast_results.df, (End_of_alignment_in_query - Start_of_alignment_in_query+1)/Query_length)*100),2)
# blast_results.df$Percent_subject_aligned <- round(abs(with(blast_results.df, (End_of_alignment_in_subject - Start_of_alignment_in_subject+1)/Subject_length)*100),2)
# head(blast_results.df)

# blast_results.df <- blast_results.df %>% dplyr::filter(Percent_query_aligned >= 50.0, Percent_subject_aligned >= 50.0)

# Query_coverage_per_HSP
# Subject_coverage_per_HSP

 # prodigal >
#   count genes >
#   BLAST genes against IMG VR >
#   filter significant hits >
#   count number of hits >
#   if >= 30% of proteins for a query have significant hit >
  # resolve consensus affilitation


