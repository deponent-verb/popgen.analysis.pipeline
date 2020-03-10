pacman::p_load(popgen.tools,tidyverse)
#Process simulation objects into dataframe

#Read all rds files in a directory

names= list.files(pattern=".rds")
genomes = lapply(names, readRDS)

#generate dataframe

#check SNP distribution ----

#The point of this section is to work out how many SNPs we should retain.
snp_dist<-snp_count(genomes)

#check snp distribution using boxplots
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

ggplot(data=snp_dist, aes(x=SNP, color=sweep_type))+ geom_density()

#bimodal distribution due to two different selection coefficients. 

temp<-snp_dist %>% filter(sweep_type=="hard") %>% filter (s==0.01) %>% select(SNP)
low_mean<-mean(temp$SNP) %>% round()
low_std<-sd(temp$SNP) %>% round()
snp_cutoff<-low_mean-2*low_std #2240

#make dataframe

df<-generate_df(sim_list = genomes,nwins = 11,snp=snp_cutoff,form="wide")
apply(df, 2, function(x) any(is.na(x)))



