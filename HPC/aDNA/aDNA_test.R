#script to make small datasets for testing purposes

.libPaths(c("/hpcfs/users/a1708050/local/RLibs",.libPaths()))

library("popgen.tools")
library("tidyverse")
library("foreach")
library("parallel")

libs = .libPaths(c("/hpcfs/users/a1708050/local/RLibs",.libPaths()))


slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (is.numeric(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
} else {
  cores = detectCores() # Figure out how many cores there are
}
cl<-makeCluster(cores)

setwd("/hpcfs/users/a1708050/mphil/ml_review/ancient_data/constant_pop")

#set DNA aging parameters
missing_rate = 0
trans_prop = 0.776
dmg_rate = c(0,0.05)
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute = c("zero","random")

#randomly split simulations into chunks for parallel SS computation
all_names = list.files(pattern=".rds")
#This line is to randomly downsample the data for testing purposes
#all_names = sample(all_names, size = 1000, replace = F)
all_names = all_names[1:1200]
n = length(all_names)/200
ID = seq(1,length(all_names),by=1)
set.seed(2)
sim_groups = split(all_names, as.factor(1:n))
ID_groups = split(ID, as.factor(1:n))

doParallel::registerDoParallel(cl,cores = cores)

#make sure you check arguments within loop
df = foreach( imp = 1:length(impute)) %:%
  foreach (r = 1:length(dmg_rate)) %:%
  foreach(i = 1:length(sim_groups)) %dopar% {
    
    #ensure correct library and directory for each core
    .libPaths(libs)
    setwd("/hpcfs/users/a1708050/mphil/ml_review/ancient_data/constant_pop")
    
    #load a small set of 100 simulations
    genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
    genomes = unlist(genomes, recursive = F)
    
    #compute SS on the small set
    popgen.tools::ancient_generate_df(sim_list = genomes,nwins = 5,
                                      split_type="base",trim_sim = F,missing_rate = 0,
                                      trans_prop = trans_prop,dmg_rate = dmg_rate[r],ascertain_indices = asc_indices,
                                      impute_method = impute[imp],ID = ID_groups[[i]])
    
    #remove the simulations from memory once we finished computing SS
    #rm(genomes)
  }

#need to unlist each of the foreach loops
df = unlist(df, recursive = F)
df = unlist(df, recursive = F)

final_df = data.table::rbindlist(df, use.names = T, fill = F, idcol = T)
readr::write_csv(final_df, file ="/hpcfs/users/a1708050/mphil/ml_review/ancient_data/dataframes/test_set2.csv")

