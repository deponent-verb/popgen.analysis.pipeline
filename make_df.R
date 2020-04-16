pacman::p_load(popgen.tools,tidyverse)
#Process simulation objects into dataframe

#Read all rds files in a directory

setwd("~/work/MPhil/ml_review/data/bottleneck_sims(hard)/")
names= list.files(pattern=".rds")[1:10]
genomes = lapply(names, readRDS)

#generate dataframe

#check SNP distribution ----

#The point of this section is to work out how many SNPs we should retain.
snp_dist<-snp_count(genomes)
anyNA(snp_dist)

#check snp distribution using boxplots
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

ggplot(data=snp_dist, aes(x=SNP, color=sweep_type))+ geom_density()

#bimodal distribution due to two different selection coefficients. 

temp<-snp_dist %>% filter(sweep_type=="hard") %>% filter (s==0.1) %>% select(SNP)
low_mean<-mean(temp$SNP) %>% round()
low_std<-sd(temp$SNP) %>% round()
snp_cutoff<-low_mean-2*low_std #2240

#let's just try 1000 for now
snp_cutoff=1000

####
doParallel::registerDoParallel()
a=Sys.time()
df<-generate_df(sim_list = genomes,nwins = 11,
                split_type="base",snp=1000,form="wide")
b=Sys.time()
####
#make dataframe

doParallel::registerDoParallel()
a=Sys.time()
df<-generate_df(sim_list = genomes,nwins = 11,snp=snp_cutoff,form="wide",
                LD_downsample = T, ds_prop = 0.1, ds_seed = 1)
b=Sys.time()
apply(df, 2, function(x) any(is.na(x)))
#readr::write_csv(df,path="~/Documents/GitHub/popgen.analysis.pipeline/data/bottleneck.csv")
readr::write_csv(df,path="~/work/MPhil/data/dataframes/constant_pop.csv")

#combining constant pop with bottleneck simulations
cpop<-read_csv("~/work/MPhil/data/dataframes/constant_pop.csv")
bt<-read_csv("~/work/MPhil/data/dataframes/only_bottleneck.csv")
df<-rbind(cpop,bt)
readr::write_csv(df,path="~/Documents/GitHub/popgen.analysis.pipeline/data/bt_cpop.csv")
