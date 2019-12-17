#load packages
pacman::p_load("popgen.tools","tidyverse","ggplot2","GGally")

#read in data
#hard<-readRDS("~/work/MPhil/data/hard.rds")
#soft<-readRDS("~/work/MPhil/data/soft.rds")
#df<-c(hard,soft)
df<-readRDS("~/work/MPhil/data/toy_data.rds")

#snp distribution----

#check SNP distribution
#snp_dist<-bind_rows(snp_count(hard),snp_count(soft))
snp_dist<-snp_count(df)


#convert rows with s=0 to neutral
snp_dist<-snp_dist %>% mutate(sweep_type=ifelse(s==0,"neutral",sweep_type))
snp_dist$sweep_type<-snp_dist$sweep_type %>% as.factor()

#check snp distribution
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

snp_dist %>% filter(sweep_type=="1") %>% summary()

#hs_dist=snp_dist %>% filter(sweep_type=="hard") %>% select(SNP)
ggplot(data=snp_dist, aes(x=SNP, color=sweep_type))+ geom_density()

#bimodal distribution due to two different selection coefficients. 

temp<-snp_dist %>% filter(sweep_type=="1") %>% filter (s==0.01) %>% select(SNP)
low_mean<-mean(temp$SNP) %>% round()
low_std<-sd(temp$SNP) %>% round()
snp_cutoff<-low_mean-2*low_std #1357

#generate the dataframe ----

#data<-c(hard,soft)
#saveRDS(data,file="~/work/MPhil/data/toy_data.rds")

data<-readRDS("~/work/MPhil/data/toy_data.rds")

df<-generate_df(sim_list = data,win_split = 10,snp=snp_cutoff)
write_csv(df,path="./data/toy_df.csv")

## Read in dataframe with raw data

df<-read_csv("./data/toy_df.csv")

df<-as_tibble(df)
df$sweep<-df$sweep %>% as.factor()
str(df)

df<-df %>% drop_na()

#sanity checking of summary statistics (raw data) ----

#index of starting column with summary stats
sum_start=4

#number of splits in a window
wins=11

ggparcoord(data=df, columns=15:25,groupColumn = 1)

#unscaled parcoords
for(i in 1:5){
  p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
  print(p)
}

i=1
 p<-ggparcoord(data=df%>%filter(sweep=="neutral"),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df%>%filter(sweep=="soft"),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df%>%filter(sweep=="hard")%>%filter(s_coef==0.01),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 temp<-df%>%filter(sweep=="hard") %>% head(50)
 p<-ggparcoord(data=temp,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 plot(p)

#scaled parcoords
for(i in 1:5){
  p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1,scale="center")
  print(p)
}

#check boxplots

ggplot(data=df,aes(x =sweep,y=D_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h1_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h2_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h12_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h123_6))+geom_boxplot()+scale_y_log10()



#### outdated code
i<-rep(0:5)
ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)

df %>% pivot_longer(H1:h123_10) %>% 
  mutate(batch = str_remove_all(name, '\\d')) %>% 
  ggplot(aes(name, value)) + geom_line(aes(group = ID)) + 
  facet_wrap(~batch, scales = "free")

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

