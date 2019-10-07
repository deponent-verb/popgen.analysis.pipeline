pacman::p_load("popgen.tools")

data<-readRDS("./data/toy_set.rds")

df<-generate_df(sim_list = data,win_split = 10)

sim<-data[[414]]
sub_win(sim$genomes,2)

#Reading in large dataset ----
hard<-readRDS("~/work/MPhil/data/hard.rds")
soft<-readRDS("~/work/MPhil/data/soft.rds")
neutral<-readRDS("~/work/MPhil/data/neutral.rds")

data<-list(hard,neutral,soft)
data<-unlist(data,recursive = F)
saveRDS(data,"~/work/MPhil/data/fiveK.rds")

df<-readRDS("~/work/MPhil/data/fiveK.rds")

df<-generate_df(sim_list = df,win_split = 10)

