#script to check SFS of ancient simulations

library("popgen.tools")

setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")

#set DNA aging parameters
missing_rate = 0.1
trans_prop = 0.776
dmg_rate = 0.05
asc_indices = lapply( seq(99,119,by=2), function(d){c(d,d+1)})
impute_method = "zero"

#for testing purposes on home machine
setwd("~/work/MPhil/ml_review/ancient_data/constant_pop/")


#randomly split simulations into chunks for parallel SS computation
all_names = list.files(pattern=".rds")
#This line is to randomly downsample the data for testing purposes
#all_names = sample(all_names, size = 50, replace = F)
#all_names = all_names[1:10]
#n = length(all_names)/300
# ID = seq(1,length(all_names),by=1)
# set.seed(2)
# sim_groups = split(all_names, as.factor(1:n))
# ID_groups = split(ID, as.factor(1:n))

sfs = list()

for(i in all_names){
  print(i)
  genome = readRDS(i)
  sfs[[i]] = compute_sfs(sim = genome,nwins = 5,missing_rate = missing_rate,trans_prop = trans_prop,
              dmg_rate = dmg_rate,ascertain_indices = asc_indices,impute_method = impute_method)
}

saveRDS(sfs, file = "~/Desktop/sfs.rds")



# genomes = lapply(sim_groups[[i]], function(d){ lapply(d,readRDS)}) 
# temp <- lapply(all_names, function(S){compute_sfs(sim = S,nwins = 5,missing_rate = missing_rate,trans_prop = trans_prop,
#                                           dmg_rate = dmg_rate,ascertain_indices = asc_indices,impute_method = impute_method)} )
