#common data cleaning script to be used across all workflows (for modern data)
library(tidyverse)

#read in data
cpop_genomes<-read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/split_snp/snp_cpop.csv")
btl_genomes <- read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/split_snp/snp_btl.csv")
snp_set = bind_rows(cpop_genomes,btl_genomes)

#X1 is the actual ID label
snp_set <- subset(snp_set, select = -c(ID, .id))
snp_set <- tibble::rowid_to_column(snp_set, "ID")

#create a new factor variable to represent the different demographic scenarios

snp_set <- snp_set %>%
  mutate(demography = paste0("t1:",bottle_time1,"_s:",bottle_size1,"_t2:",
                             bottle_time2)) %>%
  mutate(severity = (bottle_time1-bottle_time2)*bottle_size1)

snp_set <- subset(snp_set,
                  select = -c(bottle_time1:bottle_size2))

#convert the null bottleneck into cpop for clarity and ease of reading
null_index <- snp_set$demography == "t1:0_s:1_t2:0"
snp_set$demography[null_index] <- 'cpop'
snp_set$demography <- as.factor(snp_set$demography)

#remove LD_max and LD_avg (made up bad stats)

snp_set = subset(snp_set,
                 select = -c(LD_avg_1:LD_max_11))

#remove block information
snp_set = subset(snp_set,
                 select = -c(block_base_length_1:block_snp_length_11))

write_csv(snp_set,"~/work/MPhil/ml_review/data/hubs_data/dataframes/0.25ds_base_set.csv")
write_csv(snp_set,"~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_chunky1_data.csv")
