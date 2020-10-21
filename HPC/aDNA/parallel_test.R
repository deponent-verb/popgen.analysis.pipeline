.libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))

library("popgen.tools")
library("tidyverse")
library("foreach")
library("parallel")
library("doSNOW")

libs = .libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))


slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (is.numeric(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
} else {
  cores = detectCores() # Figure out how many cores there are
}
cl<- makeCluster(cores)

setwd("/fast/users/a1708050/mphil/ml_review/ancient_data/constant_pop")

#set DNA aging parameters
missing_rate = seq(0,0.95,by=0.1)
#missing_rate = 0
trans_prop = 0.776
#dmg_rate = seq(0.01,0.05,by=0.01)
dmg_rate = 0.05
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute = c("zero","random")

#for testing purposes on home machine
#setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")


#randomly split simulations into chunks for parallel SS computation
all_names = list.files(pattern=".rds")
set.seed(2)
#This line is to randomly downsample the data for testing purposes
#all_names = sample(all_names, size = 1000, replace = F)
n = length(all_names)/500
sim_groups = split(all_names, as.factor(1:n))

doParallel::registerDoParallel(cl,cores = cores)

df = foreach (r = 1:length(missing_rate)) %:%
  foreach( imp = 1:length(impute)) %:%
  #df = foreach (r = 1:length(dmg_rate)) %:%
  foreach(i = 1:length(sim_groups)) %dopar% {
    
    print(impute)
    #ensure correct library and directory for each core
    #.libPaths(libs)
    #setwd("/fast/users/a1708050/mphil/ml_review/ancient_data/constant_pop")
    
    #load a small set of 100 simulations
    # genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
    # genomes = unlist(genomes, recursive = F)
    # 
    # #compute SS on the small set
    # popgen.tools::ancient_generate_df(sim_list = genomes,nwins = 5,
    #                                   split_type="mut",trim_sim = F,missing_rate = missing_rate[r],
    #                                   trans_prop = trans_prop,dmg_rate = dmg_rate,ascertain_indices = asc_indices,
    #                                   impute_method = impute[imp])
    
    #remove the simulations from memory once we finished computing SS
    #rm(genomes)
  }
