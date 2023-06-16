

###########################################################################
###########################################################################
#        NOTE -- THIS IS NOT NEEDED WHEN USING DATA FROM ENA
###########################################################################
###########################################################################

# 0. Prep environment 
# 1. Remove Primers

#########################
# 0. Prep environment 
#########################

#load library
library(dada2); packageVersion("dada2")
library(tidyverse)
cutadapt <- "/Users/fmazel/Softwares/ls/envs/cutdapt_env/bin/cutadapt" # path to cutadapt on the computer 


# location of data on computer 
path <- "/Users/fmazel/Data/Criquets/p419_run191016_16S_RawData_Results_Reports/p419_run191016_16S_raw_data/a_data/gz" # CHANGE ME to the directory containing the fastq files after unzipping.

# location of F and R reads on computer 
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# summary of pathways 
reads_path = tibble(ForwardReadsPath = fnFs,
                    ReverseReadsPath = fnRs,
                    ForwardReadsName = sort(list.files(path, pattern="_R1_001.fastq", full.names = F)),
                    ReverseReadsName = sort(list.files(path, pattern="_R2_001.fastq", full.names = F)))

path_to_trimmed_reads = "/Users/fmazel/Data/Criquets/MyWork/TrimmedRawReads_error_Submitted/"

reads_path = reads_path %>% 
  mutate(Trimmed_ForwardReadsPath = paste0(path_to_trimmed_reads,ForwardReadsName),
         Trimmed_ReverseReadsPath = paste0(path_to_trimmed_reads,ReverseReadsName))


#########################
# 1. Remove Primers
#########################

system2(cutadapt, args = "--help")
system2(cutadapt, args = "--version") # 4.4 
R_primers_Freads = "GACTACHVGGGTATCTAATCC"
F_primers_Rreads = "CCTACGGGNGGCWGCAG"
error=0.15
error*nchar(R_primers_Freads) # 3.15 errors 
error*nchar(F_primers_Rreads) # 2.55 errors

for (i in 1:404)
{
  print(i)
  params_cutadapt <- c('-e', error, # 
                       '-m', 100, # discard all reads < 100 pb
                       '-g', R_primers_Freads, # Forward barcode 
                       '-G', F_primers_Rreads,  # Reverse barcode 
                       '-o',reads_path$Trimmed_ForwardReadsPath[i], # R1 output
                       '-p',reads_path$Trimmed_ReverseReadsPath[i], # R2 output 
                       reads_path$ForwardReadsPath[i], # Input R1
                       reads_path$ReverseReadsPath[i] # Input R2 
  )
  
  
  system2(cutadapt, args = params_cutadapt)
  
}

