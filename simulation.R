pacman::p_load(popgen.tools)
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

bottlenecks = tibble::tibble(size=c(0.5,1),time=c(5,10))

# selection coefficients ----

selection=c(0,0.001, 0.01)

## Simulation Loops ----

nsim=1000
setwd("~/work/MPhil/data/batch_data/")
sweep_type="hard"

#loop for hard sweeps

doParallel::registerDoParallel()
a=Sys.time()
for(s in selection){
#  print(s)
  for(i in 1:nsim){
#    print(s)
    sim = discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                genome_length = nBases, samplesize = samplesize,
                s = s, discoal_path = discoal_path,
                sweep=sweep_type, fix_time = 1)
    name = paste("hardsim_s",s,"_n",i,".rds",sep="")
#    print(name)
    saveRDS(sim,file=name)
  }
}
b=Sys.time()

#1000 sims each, took ~3.6 hours




