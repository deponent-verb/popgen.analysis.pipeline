#PCA workflow

#This script goes through each demographic model and fits a PCA. 
#We will use this to make some plots. In particular, we want to 
#investigate which summary stats are differentiating the sweeps
#for the different demographies. 

library(tidymodels)
library(tidyverse)

#load data
genomes <- read_csv("./data/bt_cpop.csv")
genome_data <- genomes %>%
  select(sweep, H_1:h123_11, demography)

pca_rec <- recipe(data = place_holder, sweep ~.) %>% #set our group to be sweep
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep() 


pca_wkfl <- workflow()