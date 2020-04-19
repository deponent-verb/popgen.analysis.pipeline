pacman::p_load(popgen.tools,parallel,foreach)
#simulate training data for analysis

# Species Model (Human) ----

mu=1.5e-8
recomb_rate=1e-8
Ne=1e4
nBases=1e6
samplesize=100
discoal_path="~/work/programs/discoal/discoal"

#sweeps are recently fixed (1 generation from sampling)
fix=1
#starting frequency for soft sweeps 
start_f=0.1

# Demographic Models ----

# Constant popsize. Leave out popsize_changes param. 

# Bottleneck scenarios.

bottlenecks = list()
recovery = c(1600, 8000)
duration = c(80, 800, 8000)
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

nsim=600
setwd("~/work/MPhil/data/bottleneck_sims(neutral)/")
sweep_type="hard"

#loop for hard sweeps, with bottlenecks

doParallel::registerDoParallel()
a=Sys.time()
for(s in selection){
  #  print(s)
  for(b in bottlenecks){
    for(i in 1:nsim){
      # print(b$size[1])
      # print(b$time[1])
      # print(b$time[2])
      # print(s)
      sim = discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                        genome_length = nBases, samplesize = samplesize,
                        s = s, discoal_path = discoal_path,
                        sweep=sweep_type, fix_time = 1, popsize_changes = b)
      name = paste("hardsim_s",s,"_n",i,"_b",b$size[1],"_t1",b$time[1],"_t2",b$time[2],
                   ".rds",sep="")
      print(name)
      saveRDS(sim,file=name)
    }
  }
}
b=Sys.time()

# constant popsize simulations

nsim=600
setwd("~/work/MPhil/ml_review/data/constantpop/")
sweep_type="hard"

cores=detectCores()
cl<-makeCluster(cores)

doParallel::registerDoParallel(cl,cores = cores)

a=Sys.time()
foreach(s = 1:length(selection)) %dopar% {
    for(i in 101:(nsim+100)) {
      sim = popgen.tools::discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                        genome_length = nBases, samplesize = samplesize,
                        s = selection[s], discoal_path = discoal_path,
                        sweep=sweep_type, fix_time = 1)
      name = paste("hardsim_s",selection[s],"_n",i,"constant_pop",
                   ".rds",sep="")
      print(name)
      saveRDS(sim,file=name)
    }
}

b=Sys.time()

# foreach(s = 1:length(selection)) %dopar% {
#   print(selection[s])
# }



