pacman::p_load(popgen.tools,parallel,foreach)
#simulate training data for analysis

# Species Model (Human) ----

mu=1.5e-8
recomb_rate=1e-8
Ne=1e4
nBases=1e6
samplesize=100
discoal_path="~/work/programs/discoal/discoal"

#sweeps are recently fixed at time of sampling
fix=0
#starting frequency for soft sweeps 
start_f=0.1

# Demographic Models ----

# Constant popsize. Leave out popsize_changes param. 

# Bottleneck scenarios.

bottlenecks = list()
duration = c(1600, 8000) #t1 = duration + t2
recovery = c(80, 800, 8000) #t2
strength = c(0.05, 0.1, 0.5)

counter=1
for(r in recovery){
  for(d in duration){
    for(s in strength){
      bottlenecks[[counter]] = tibble::tibble(size=c(s,1),time=c(d+r,r))
      counter=counter+1
    }
  }
}

# selection coefficients ----

selection=c(0,100,250,500,750,1000,2000)/(2*Ne)
#selection=0

## Simulation Loops ----

nsim=1000
setwd("~/work/MPhil/ml_review/data/hubs_data/bottlenecks/")
sweep_type="hard"

#loop for hard sweeps, with bottlenecks

cores=detectCores()
cl<-makeCluster(cores)
doParallel::registerDoParallel(cl,cores = cores)

a=Sys.time()
foreach(s = 1:length(selection)) %dopar%{
  for(b in bottlenecks){
    for(i in 1:nsim){
      sim = popgen.tools::discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                        genome_length = nBases, samplesize = samplesize,
                        s = selection[s], discoal_path = discoal_path,
                        sweep=sweep_type, fix_time = fix, popsize_changes = b)
      name = paste("hardsim_s",selection[s],"_n",i,"_b",b$size[1],"_t1",b$time[1],"_t2",b$time[2],
                   ".rds",sep="")
      print(name)
      saveRDS(sim,file=name)
    }
  }
}
b=Sys.time()

#1000 sims for 18 bottlenecks, 7 s_coef took 1.33 days on 4 cores

# constant popsize simulations ----

nsim=1000
setwd("~/work/MPhil/ml_review/data/hubs_data/constant_pop/")
sweep_type="hard"

cores=detectCores()
cl<-makeCluster(cores)

doParallel::registerDoParallel(cl,cores = cores)

a=Sys.time()
foreach(s = 1:length(selection)) %dopar% {
    for(i in 1:nsim) {
      sim = popgen.tools::discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                        genome_length = nBases, samplesize = samplesize,
                        s = selection[s], discoal_path = discoal_path,
                        sweep=sweep_type, fix_time = fix)
      name = paste("hardsim_s",selection[s],"_n",i,"constant_pop",
                   ".rds",sep="")
      print(name)
      saveRDS(sim,file=name)
    }
}
b=Sys.time()

#20mins to do 100, for each 7 selection coef
