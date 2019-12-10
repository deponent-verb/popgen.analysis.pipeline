#load packages
pacman::p_load("popgen.tools","tidyverse","ggplot2","GGally")

#read in data
hard<-readRDS("~/work/MPhil/data/hard.rds")
soft<-readRDS("~/work/MPhil/data/soft.rds")

#snp distribution----

#check SNP distribution
snp_dist<-bind_rows(snp_count(hard),snp_count(soft))

#convert rows with s=0 to neutral
snp_dist<-snp_dist %>% mutate(sweep_type=ifelse(s==0,"neutral",sweep_type))
snp_dist$sweep_type<-snp_dist$sweep_type %>% as.factor()

#check snp distribution
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

ggplot(snp_dist,aes(sweep_type=="hard",SNP))+geom_boxplot()

snp_dist %>% filter(sweep_type=="hard") %>% summary()

#hs_dist=snp_dist %>% filter(sweep_type=="hard") %>% select(SNP)
ggplot(data=snp_dist, aes(x=SNP, color=sweep_type))+ geom_histogram(bins=50)

#discussion. Areas of each group don't look the same. Bimodal distribution for hard sweeps. 



#snp_include<-temp$SNP %>% mean() %>% round()
#handwave 1600, change later
snp_include=1600

#generate the dataframe ----

#data<-c(hard,soft)
#saveRDS(data,file="~/work/MPhil/data/toy_data.rds")

data<-readRDS("~/work/MPhil/data/toy_data.rds")

df<-generate_df(sim_list = data,win_split = 10,snp=snp_include)
write_csv(df,path="./data/toy_df.csv")

df<-read.csv("./data/toy_df.csv")
df<-as_tibble(df)
str(df)

#sanity checking of summary statistics

#index of starting column with summary stats
sum_start=4

#number of splits in a window
wins=10

i<-rep(0:5)
f<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+i*wins-1),groupColumn=1)



ggplot(data=df, aes(x=dist,y=D)) + geom_line(aes(color=sweep)) + scale_y_continuous(trans="log10")


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

