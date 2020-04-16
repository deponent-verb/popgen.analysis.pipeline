.libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))

library("popgen.tools")
library("tidyverse")
library("parallel")

slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (is.numeric(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
} else {
  cores = detectCores() # Figure out how many cores there are
}
cl<-makeCluster(cores)

#Generates dataframe out of simulation objects

#Read all rds files in a directory

setwd("/fast/users/a1708050/mphil/ml_review/data/bottleneck_sims(hard)")
names= list.files(pattern=".rds")[1:3]
genomes = lapply(names, readRDS)

doParallel::registerDoParallel(cl,cores = 8)
a=Sys.time()
df<-generate_df(sim_list = genomes,nwins = 11,
                split_type="base",snp=1000,form="wide")
b=Sys.time()
readr::write_csv(df,path="/fast/users/a1708050/mphil/ml_review/data/dataframes/sbottle_necks.csv")
b-a
