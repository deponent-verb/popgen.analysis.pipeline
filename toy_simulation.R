pacman::p_load(popgen.tools,parallel)

#population params----
mu=1.5e-8
recomb_rate=1e-8
Ne=1e4
nBases=1e6
samplesize=20
s=c(0,0.001,0.01)
#s=c(10,50,100,500,1000)*(1/(2*Ne))
fix=1
discoal_path="~/work/programs/discoal/discoal"

#####################################
#batch simulation function----
#takes a selection coefficients and runs n hard sweep simulations

#input: N the number of simulations per selection coefficient
#s: a selection coefficient

#return: a list of sim_objs
batch_sim<-function(select_coeff,N,sweep_type){
  sims<-list()
  for(i in 1:N){
    sims[[i]]<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=select_coeff,discoal_path=discoal_path,fix_generation=fix,sweep=sweep_type)
  }
  return(sims)
}

#Running simulations----

#mclapply takes first element, runs it on a core. And so on. 
num_sim=1000

Sys.time()
hard=mclapply(s,batch_sim,N=num_sim,sweep_type="hard",mc.cores=4)
Sys.time()
hard<-unlist(hard,recursive = F)
saveRDS(hard,"~/work/MPhil/data/hard.rds")
#8mins for num_sim=100


Sys.time()
soft=mclapply(s,batch_sim,N=num_sim,sweep_type="soft",mc.cores=4)
Sys.time()
soft<-unlist(soft,recursive = F)
saveRDS(soft,"~/work/MPhil/data/soft.rds")
#1min for num_sim 100

# Sys.time()
# neutral<-mclapply(s,batch_sim,N=num_sim,sweep_type="neutral",mc.cores=4)
# Sys.time()
# neutral<-unlist(neutral,recursive = F)
# saveRDS(neutral,"~/work/MPhil/data/neutral.rds")
#10mins for num_sim 100

#debugging
#neutral<-lapply(s,batch_sim,N=num_sim,sweep_type="neutral")


##### Putting everything together into one R object list
# df<-list(hard,soft,neutral)
# df<-unlist(df,recursive = F)
# saveRDS(df,"~/work/MPhil/data/toy_set.rds")

##load data
df<-readRDS("~/work/MPhil/data/toy_data.rds")
data<-generate_df(sim_list = df,win_split = 11,snp=1357)

