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
cl<-makeCluster(cores)

cores

a=Sys.time()
#Generates dataframe out of simulation objects

#Read all rds files in a directory

setwd("/fast/users/a1708050/mphil/ml_review/hubsdata/bottlenecks")
i=3
all_names = list.files(pattern=".rds")
all_names = all_names[(42000*(i-1)+1) :(42000*i)]

set.seed(1)
n = length(all_names)/500
sim_groups = split(all_names, as.factor(1:n))


doParallel::registerDoParallel(cl,cores = cores)

df = foreach(i = 1:length(sim_groups)) %dopar% {
  
  #ensure correct library and directory for each core
  .libPaths(libs)
  setwd("/fast/users/a1708050/mphil/ml_review/hubsdata/bottlenecks")
  
  #load a small set of 100 simulations
  genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
  genomes = unlist(genomes, recursive = F)
  
  #compute SS on the small set
  popgen.tools::generate_df(sim_list = genomes,nwins = 11,
                            split_type="base",trim_sim = F,form="wide",
                            LD_downsample = T, ds_prop = 0.25)
  
  #remove the simulations from memory once we finished computing SS
  #rm(genomes)
}

final_df = data.table::rbindlist(df, use.names = T, fill = F, idcol = T)
readr::write_csv(final_df, path= "/fast/users/a1708050/mphil/ml_review/hubsdata/dataframes/notrim_base_btl3.csv")