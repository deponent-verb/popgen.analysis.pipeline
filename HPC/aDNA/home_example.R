library("popgen.tools")
library("tidyverse")
library("foreach")
library("parallel")
#library("doSNOW")

setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")

#set DNA aging parameters
missing_rate = seq(0,0.95,by=0.1)
#missing_rate = 0
trans_prop = 0.776
#dmg_rate = seq(0.01,0.05,by=0.01)
dmg_rate = 0.05
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute = c("zero","random")
denoise_method = c("none","majority_flip","cluster")


#for testing purposes on home machine
setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")


#randomly split simulations into chunks for parallel SS computation
all_names = list.files(pattern=".rds")
#This line is to randomly downsample the data for testing purposes
#all_names = sample(all_names, size = 50, replace = F)
all_names = all_names[1:1200]
n = length(all_names)/300
ID = seq(1,length(all_names),by=1)
set.seed(2)
sim_groups = split(all_names, as.factor(1:n))
ID_groups = split(ID, as.factor(1:n))

#setup strategy stops mac from hanging
cl<-parallel::makeCluster(4,setup_strategy = "sequential")
doParallel::registerDoParallel(cl,cores = cores)

start = Sys.time()
df = foreach (r = 1:length(missing_rate)) %:%
  foreach( imp = 1:length(impute)) %:%
  foreach (d = 1:length(denoise_method)) %:%
  foreach(i = 1:length(sim_groups)) %dopar% {
    
    #ensure correct library and directory for each core
    setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")
    
    #load a small set of 100 simulations
    genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
    genomes = unlist(genomes, recursive = F)
    
    #compute SS on the small set
    popgen.tools::ancient_generate_df(sim_list = genomes,nwins = 5,
                                      split_type="mut",trim_sim = F,missing_rate = missing_rate[r],
                                      trans_prop = trans_prop,dmg_rate = dmg_rate,ascertain_indices = asc_indices,
                                      impute_method = impute[imp],ID = ID_groups[[i]],
                                      denoise_method = denoise_method[d])
    
    #remove the simulations from memory once we finished computing SS
    #rm(genomes)
  }
#6:27
end = Sys.time()
#need to unlist each of the foreach loops
# df = unlist(df, recursive = F)
# df = unlist(df, recursive = F)
temp = df

temp = unlist(temp, recursive = F)
temp = unlist(temp, recursive = F)
temp = unlist(temp, recursive = F)

final_df = data.table::rbindlist(temp, use.names = T, fill = F, idcol = T)
readr::write_csv(final_df,path="~/work/MPhil/ml_review/ancient_data/ancient_df/aDNA_small.csv")
