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
missing_rate = c(0,0.2,0.4,0.6,0.8,0.95)
trans_prop = 0.776
dmg_rate = 0.05
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute = c("zero","random")
denoise = c("none", "cluster")
sim_batch = 100

#for testing purposes on home machine
setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")


#randomly split simulations into chunks for parallel SS computation
all_names = list.files(pattern=".rds")
#shuffle the simulations
set.seed(3)
all_names = sample(all_names, size = length(all_names), replace = F)
#This line is to randomly downsample the data for testing purposes
# set.seed(22)
#=all_names = sample(all_names, size = 1000, replace = F)
n = length(all_names)/sim_batch
ID = seq(1,length(all_names),by=1)
set.seed(2)
sim_groups = split(all_names, as.factor(1:n))
set.seed(2)
ID_groups = split(ID, as.factor(1:n))

doParallel::registerDoParallel(cl,cores = cores)

#home test
#cl<-parallel::makeCluster(4,setup_strategy = "sequential")

df = foreach( imp = 1:length(impute)) %:%
  foreach (d = 1:length(denoise)) %:%
  foreach(i = 1:length(sim_groups)) %dopar% {
    
    #ensure correct library and directory for each core
    .libPaths(libs)
    setwd("/hpcfs/users/a1708050/mphil/ml_review/ancient_data/constant_pop")
    
    #test for home machine. Comment out for Phoenix. 
    # library(popgen.tools)
    # setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")
    
    
    #load a small set of 100 simulations
    genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
    genomes = unlist(genomes, recursive = F)
    
    #randomly sample one missing rate
    miss = sample(missing_rate, size = 1)
    
    #compute SS on the small set
    popgen.tools::ancient_generate_df(sim_list = genomes,nwins = 5,
                                      split_type="mut",trim_sim = F,missing_rate = miss,
                                      trans_prop = trans_prop,dmg_rate = dmg_rate,ascertain_indices = asc_indices,
                                      impute_method = impute[imp],ID = ID_groups[[i]],denoise_method = denoise[d])
    
    #remove the simulations from memory once we finished computing SS
    #this line needs to be commented out or we get null dataframes
    #rm(genomes)
  }

#need to unlist each of the foreach loops
df = unlist(df, recursive = F)
df = unlist(df, recursive = F)

final_df = data.table::rbindlist(df, use.names = T, fill = F, idcol = T)
head(final_df)
readr::write_csv(final_df, file ="/hpcfs/users/a1708050/mphil/ml_review/ancient_data/dataframes/aDNA_cpop.csv")