#script to make small datasets for testing purposes

.libPaths(c("/hpcfs/users/a1708050/local/RLibs",.libPaths()))

library("parallel")

libs = .libPaths(c("/hpcfs/users/a1708050/local/RLibs",.libPaths()))


slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (is.numeric(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
} else {
  cores = detectCores() # Figure out how many cores there are
}
cl<-makeCluster(cores)