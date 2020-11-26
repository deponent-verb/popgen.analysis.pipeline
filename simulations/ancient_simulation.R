pacman::p_load(popgen.tools,parallel,foreach)

#simulate training data for analysis

# Species Model (Human) ----

mu=1.5e-8
recomb_rate=1e-8
Ne=1e4
nBases=1e6
samplesize=120
discoal_path="~/work/programs/discoal/discoal"
fix=0 #sweeps are recently fixed at time of sampling
selection=c(0,100,250,500,750,1000,2000)/(2*Ne)

#population tree
demes = 2
sample_dist = c(100,20)
deme_join = tibble::tibble(time = 50000/25,pop1 = 0, pop2=1)

nsim=500*6
setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")
sweep_type="hard"

cores=detectCores()
cl<-makeCluster(cores,setup_strategy = "sequential")

doParallel::registerDoParallel(cl,cores = cores)

a=Sys.time()
foreach(s = 1:length(selection)) %dopar% {
  for(i in 1:nsim) {
    sim = popgen.tools::discoal_sim(mu=mu,recomb_rate = recomb_rate, Ne = Ne,
                                    genome_length = nBases, samplesize = samplesize,
                                    s = selection[s], discoal_path = discoal_path,
                                    sweep=sweep_type, fix_time = fix,
                                    demes = demes, sample_dist = sample_dist, deme_join = deme_join)
    name = paste("ancient_hardsim_s",selection[s],"_n",i,"constant_pop",
                 ".rds",sep="")
    print(name)
    saveRDS(sim,file=name)
  }
}
b=Sys.time()

#nsim = 10, takes 2.6mins
#3000 sims took 13 hours
