#load packages
pacman::p_load("popgen.tools","tidyverse","ggplot2")

#read in data
hard<-readRDS("~/work/MPhil/data/hard.rds")
soft<-readRDS("~/work/MPhil/data/soft.rds")

#check SNP distribution
snp_dist<-bind_rows(snp_count(hard),snp_count(soft))

#convert rows with s=0 to neutral
snp_dist<-snp_dist %>% mutate(sweep_type=ifelse(s==0,"neutral",sweep_type))
snp_dist$sweep_type<-snp_dist$sweep_type %>% as.factor()

#check snp distribution
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

ggplot(snp_dist,aes(sweep_type=="hard",SNP))+geom_boxplot()

snp_dist %>% filter(sweep_type=="hard") %>% summary()

#df<-generate_df(sim_list = data,win_split = 10)



DF<-readRDS("~/work/MPhil/data/df.rds")

sim<-data[[414]]
sub_win(sim$genomes,2)


#generate new IDs so that they are not informative. 

size<-nrow(df)/10
x<-rep(1:size)
ID<-sample(x,size=size)

foo<-function(y,it){
  ans<-rep(y,it)
  return(ans)
}


new_ID<-map2(ID,10,foo)
test<-unlist(new_ID)

#The interpretibility of the model is helpful for understanding the underlying biology that's giving you that answer.
#Help work out what signatures are important for finding sweeps in genome. 

