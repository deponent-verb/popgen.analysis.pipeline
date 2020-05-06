library(profvis)
library(popgen.tools)

mu=1.5e-8
recomb_rate=1e-8
Ne=1e4
nBases=1e6
samplesize=100
s=0
fix=1
discoal_path="~/work/programs/discoal/discoal"
sweep_type="hard"
nwins=11
id=1
seed=c(1,1)

temp<-discoal_sim(mu=mu,recomb_rate=recomb_rate,Ne=Ne,genome_length=nBases,samplesize=samplesize,s=s,discoal_path=discoal_path,fix_time=fix,sweep=sweep_type,seed=seed)

#no downsampling
profvis({
  doParallel::registerDoParallel()
  output<-sum_stats(sim=temp,split_type="base",nwins = nwins, ID=id, snp=1000)
})

profvis({
  doParallel::registerDoParallel()
  output<-sum_stats(sim=temp,split_type="base",nwins = nwins, ID=id, snp=1000,
                    LD_downsample = T, ds_prop = 0.1, ds_seed = 5)
})

profvis({
  temp = foreach::foreach (i = 1:length(tune_blocks), .packages = c('tune','yardstick')) %dopar% {
    
    tune::tune_grid(meta_workflow,
                    resamples = cv_splits,
                    grid = tune_blocks[[i]], 
                    metrics=metric_set(accuracy),
                    control=control_grid(save_pred = TRUE)) 
  }
})

profvis({
  tuning = tune_grid(meta_workflow,
                     resamples = cv_splits,
                     grid = tuning_params, 
                     metrics=metric_set(accuracy),
                     control=control_grid(save_pred = TRUE))
})
