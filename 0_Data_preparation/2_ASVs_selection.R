#######################################
# Prepare dataset with dada2 pipeline #
#######################################

# 0. Prep environment 
# 2. Filter and Trim 
# 3. Error model and Run dada2
# 4. Merge reads and remove chimera
# 5. Assign TAXONOMY
# 6. Remove ASVs based on taxonomy


#########################
# 0. Prep environment 
#########################

#load library
library(dada2); packageVersion("dada2")
library(tidyverse)
library(seqinr)
library(vegan)

# location of data on computer 
path <- "/Users/fmazel/Data/Criquets/Raw_fastq_ENA/" # CHANGE ME to the directory containing the fastq files after unzipping

# path for dada2 output
path_dada2 = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/"

# ref SILVA data path 
silva_ref = "/Users/fmazel/Data/SILVA_16S/dada2_formated/silva_nr99_v138_train_set.fa.gz"

# location of F and R reads on computer 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))


#########################
# 2. Filter and Trim 
#########################

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Visual inspection of quality
plotQualityProfile(fnFs[10:12])
plotQualityProfile(fnRs[10:12])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_dada2, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_dada2, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,260),
                     maxN=0, maxEE=c(4,4), truncQ=2, rm.phix=TRUE, # maxEE=c(4,2)
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

Stats_dada2 <- tidyr::as_tibble(out,rownames="sample.names") %>% mutate(propKept=reads.out/reads.in)
summary(Stats_dada2$propKept)

plotQualityProfile(filtFs[10:12])
plotQualityProfile(filtRs[10:12])

#########################
# 3. Error model and Run
#########################

# Learn errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

saveRDS(object = errF,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/errModel_errF.Rdata")
saveRDS(object = errR,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/errModel_errR.Rdata")

# Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
head(dadaFs[[1]])

saveRDS(object = dadaFs,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/ASV_infer_dadaFs.Rdata")
saveRDS(object = dadaRs,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/ASV_infer_dadaRs.Rdata")


#########################
# 4. Merge reads and remove chimera
#########################

# Merge reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]]) # Inspect the merger data.frame from the first sample
saveRDS(object = mergers,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/mergers.Rdata")

# construct table
seqtab <- makeSequenceTable(mergers)
dim(seqtab) # 368 12127


# remove chimnera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) # 368 4915
sum(seqtab.nochim)/sum(seqtab) #  0.8942331
saveRDS(object = seqtab.nochim,file = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/seqtab.nochim.Rdata")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/track_reads_thru_dada2.txt", row.names = T)


#########################
# 5 Assign TAXONOMY
#########################

taxa138 <- assignTaxonomy(seqs=colnames(seqtab.nochim), refFasta=silva_ref, multithread=6, tryRC=T, minBoot = 60, outputBootstraps=T)
tibble_taxa138 = as_tibble(taxa138[[1]])
tibble_taxa138$seq=rownames(taxa138[[1]])
write.table(tibble_taxa138 , "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/TaxoSILVA138_RDP.txt", sep="\t", quote=F)


#########################
# 6. Remove ASVs based on taxonomy
#########################

taxonomy = tibble_taxa138 

# Remove ASVs assigned to chloroplast, mitochondria or not assigned to a Kingdom - remove blanks 
sum(taxonomy$Order=="Chloroplast", na.rm = T) #190
sum(taxonomy$Family=="Mitochondria", na.rm = T) #926
sum(is.na(taxonomy$Phylum)) #1540
sum(is.na(taxonomy$Class)) #1665
sum(is.na(taxonomy$Order)) #1665

final_taxonomy = taxonomy %>% 
  subset(!Order == "Chloroplast" | is.na(Order)) %>% 
  subset(!Family == "Mitochondria"| is.na(Family)) %>% 
  subset(!is.na(Order))


###################################
# 7. Remove ASVs based on length
###################################

final_taxonomy$length_ASV = nchar(final_taxonomy$seq) # Inspect distribution of sequence lengths
hist(final_taxonomy$length_ASV)

final_taxonomy = final_taxonomy %>% 
  subset(length_ASV>390 & length_ASV<470)
hist(final_taxonomy$length_ASV)


###################################
# 7.Choose samples 
###################################

# remove low depth samples 
seqtab.nochim_selectedASVs <- seqtab.nochim[,final_taxonomy$seq] 
dim(seqtab.nochim_selectedASVs) # 368 2025

depth <- tibble(depth=apply(seqtab.nochim_selectedASVs,1,sum),
                SampleID=rownames(seqtab.nochim_selectedASVs))

deep_enought = depth %>% 
  subset(depth>1000) %>% 
  pull(SampleID)
length(deep_enought)#336

seqtab.nochim_selectedASVs_1000reads <- seqtab.nochim_selectedASVs[deep_enought ,]
seqtab.nochim_selectedASVs_1000reads <- seqtab.nochim_selectedASVs_1000reads[,apply(seqtab.nochim_selectedASVs_1000reads,2,sum)>0] # remove ASVs without count anymore 
dim(seqtab.nochim_selectedASVs_1000reads) #  336 1957


###################################
#      9.EXPORT DATA  
###################################

# Taxonomy files 
export_taxonomy = final_taxonomy %>% 
  subset(seq %in% colnames(seqtab.nochim_selectedASVs_1000reads)) # 2053

write.table(export_taxonomy,"/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/ASV_taxonomy_filtered.txt") # 2081 ASVs
write.fasta(sequences = as.list(export_taxonomy$seq),
            names = export_taxonomy$seq,
            file.out = "/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/ASV_filtered.fasta")

# ASVs tables 
write.table(seqtab.nochim_selectedASVs_1000reads,"/Users/fmazel/Data/Criquets/Processed_data/dada2_data_submission/seqtab.nochim_selectedASVs_1000reads.txt")

