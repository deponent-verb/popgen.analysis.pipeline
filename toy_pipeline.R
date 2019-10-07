pacman::p_load("popgen.tools")

data<-readRDS("./data/toy_set.rds")

df<-generate_df(sim_list = data,win_split = 2)

sim<-data[[414]]
sub_win(sim$genomes,2)
