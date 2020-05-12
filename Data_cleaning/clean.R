#common data cleaning script to be used across all workflows

#read in data
snp_set<-read_csv("~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")

#set correct factors
snp_set$sweep = as.factor(snp_set$sweep)
snp_set$s_coef=as.factor(snp_set$s_coef)

#X1 is the actual ID label
snp_set <- subset(snp_set, select = -c(ID))

#create a new factor variable to represent the different demographic scenarios

snp_set <- snp_set %>%
  mutate(demography = paste0("t1:",bottle_time1,"_s:",bottle_size1,"_t2:",
                             bottle_time2)) %>%
  mutate(severity = (bottle_time1-bottle_time2)*bottle_size1)

#convert the null bottleneck into cpop for clarity and ease of reading
null_index <- snp_set$demography == "t1:0_s1:1_t2:0_s2:1"
snp_set$demography[null_index] <- 'cpop'
snp_set$demography <- as.factor(snp_set$demography)

write_csv(snp_set,"./data/bt_cpop.csv")
