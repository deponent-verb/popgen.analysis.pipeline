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

setwd("/fast/users/a1708050/mphil/ml_review/data/constantpop")
names = list.files(pattern=".rds")
folds = split(names, as.factor(1:cores))
genomes = lapply(folds, function(d){ lapply(d,readRDS)})


doParallel::registerDoParallel(cl,cores = cores)

df = foreach(i = 1:length(genomes)) %dopar% {
  # .libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))
  #clusterEvalQ(cl, .libPaths("/fast/users/a1708050/local/RLibs"))
  .libPaths(libs)
  
  popgen.tools::generate_df(sim_list = genomes[[i]],nwins = 11,
                            split_type="base",snp=1000,form="wide",
                            LD_downsample = T, ds_prop = 0.2)
}
b=Sys.time()

final_df = do.call(rbind,df)

readr::write_csv(final_df,path="/fast/users/a1708050/mphil/ml_review/data/dataframes/base_split/base_cpop.csv")
b-a
