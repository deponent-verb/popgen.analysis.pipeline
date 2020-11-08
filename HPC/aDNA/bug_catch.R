#bug catching script for home machine

pacman::p_load(popgen.tools)

#set DNA aging parameters
missing_rate = seq(0,0.9,by=0.1)
#missing_rate = 0
trans_prop = 0.776
#dmg_rate = seq(0.01,0.05,by=0.01)
dmg_rate = 0.05
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute = c("zero","random")
denoise_method = c("none","cluster")

setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")
set.seed(22)
all_names = list.files(pattern=".rds")
trunc_names = sample(all_names, size = 1000, replace = F)
ID = seq(1,length(trunc_names),by=1)

a = Sys.time()
df = for(r in missing_rate){
  for(imp in impute){
    for(denoise in denoise_method){
      for(i in 1:length(trunc_names)){
        print(c(r,imp,trunc_names[i],denoise))
        genome = readRDS(trunc_names[i])
        ancient_sum_stats(sim = genome, nwins = 5, split_type = "mut", 
                          trim_sim = F, missing_rate = r,ascertain_indices = asc_indices,
                          impute_method = imp,denoise_method = denoise, ID = ID[i])
      }
    }
  }
}
b = Sys.time()
#2:32

#testing individual sim objects
genome = readRDS(trunc_names[34])
ancient_sum_stats(sim = genome, nwins = 5, split_type = "mut", 
                  trim_sim = F, missing_rate = 0.9,ascertain_indices = asc_indices,
                  impute_method = "zero" ,denoise_method = "cluster", ID = 1)

# [1] "0.9"                                        
# [2] "random"                                     
# [3] "ancient_hardsim_s0.025_n601constant_pop.rds"

#goes iffy if missing rate is too high
genome = readRDS(trunc_names[2])
ancient_sum_stats(sim = genome, nwins = 5, split_type = "mut", 
                  trim_sim = F, missing_rate = 0.9,ascertain_indices = asc_indices,
                  impute_method = "random" ,denoise_method = "cluster", ID = 1)

#missing rate 0, "zero" impute, "ancient_hardsim_s0.1_n496constant_pop.rds"
#more cluster centers than data points error
genome = readRDS(trunc_names[203])
saveRDS(genome, file = "~/Desktop/test_genome.rds")
ancient_sum_stats(sim = genome, nwins = 5, split_type = "mut", 
                  trim_sim = F, missing_rate = 0,ascertain_indices = asc_indices,
                  impute_method = "zero" ,denoise_method = "cluster", ID = 1)
