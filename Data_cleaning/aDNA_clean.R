#data cleaning script for aDNA simulations

library(tidyverse)

#read in data
#ancient_genomes<-read_csv("~/work/MPhil/ml_review/ancient_data/ancient_df/ancient_cpop1.csv")
#ancient_genomes<-read_csv("~/work/MPhil/ml_review/ancient_data/ancient_df/aDNA_cpop.csv")
#ancient_genomes<-read_csv("~/work/MPhil/ml_review/ancient_data/ancient_df/aDNA_small.csv")
ancient_genomes<-read_csv("~/work/MPhil/ml_review/ancient_data/ancient_df/test_set1.csv")

#filter out unneeded params
ancient_genomes <- subset(ancient_genomes, 
                          select = -c(bottle_time1:start_freq))
ancient_genomes <- subset(ancient_genomes,
                          select = -c(block_base_length_1:block_snp_length_5))
ancient_genomes <- subset(ancient_genomes,
                          select = -c(.id))

#write_csv(ancient_genomes,"~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")

#write_csv(ancient_genomes,"~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
write_csv(ancient_genomes,"~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_testset1.csv")

