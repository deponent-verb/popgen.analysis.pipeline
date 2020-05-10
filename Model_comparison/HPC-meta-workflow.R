#Phoenix HPC version of the meta-workflow
.libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))

library(tidyverse)
library(vip)
library(tidymodels)
source("./Model_comparison/model_tune.R")
source("./Model_comparison/model_vip.R")
source("./Model_comparison/model_performance.R")

libs = .libPaths(c("/fast/users/a1708050/local/RLibs",.libPaths()))

slurm_ntasks <- as.numeric(Sys.getenv("SLURM_NTASKS")) # Obtain environment variable SLURM_NTASKS
if (is.numeric(slurm_ntasks)) {
  cores = slurm_ntasks # if slurm_ntasks is numerical, then assign it to cores
} else {
  cores = detectCores() # Figure out how many cores there are
}
cl<-makeCluster(cores)
doParallel::registerDoParallel(cl,cores = cores)

#load data from cleaning script

genomes = read_csv("./data/bt_cpop.csv")

#truncate dataset to contain response and predictors only. Used for model fitting.

genomes_SS = genomes %>% 
  dplyr::select(sweep,H_1:h123_11,demography)
skimr::skim(genomes_SS)
genomes_SS$demography <- as.factor(genomes_SS$demography)

#Partition dataset into training and testing sets.
set.seed(1066)
genome_split<-initial_split(genomes_SS,prop=0.8)
genome_train = training (genome_split)
genome_test = testing (genome_split)

#Recipe design ----

#The standard recipe just standardizes all the predictors (mean=0, var =1) 

std_recipe <- recipe(sweep ~., data=genome_train) %>% #set sweep as response variable. everything else is a predictor.
  update_role(demography, new_role = 'demography') %>% #remove demography as a predictor
  step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
  step_normalize(all_predictors()) %>% #normalize all predictors
  prep()

baked_train <- bake(std_recipe, genome_train)

#Designate model and hyperparameters ----

#Logistical regression with L2 regularisation ----
genome_lr = logistic_reg(
  mode="classification",
  penalty = tune(),
  mixture=1
) %>%
  set_engine("glmnet")

#Create set of tuning parameters
lr_grid = grid_regular(penalty(range=c(0,0.2)),
                       levels=10, 
                       original = F)

#RandomForest classifier

genome_rf<-rand_forest(
  mode="classification",
  mtry=tune(),
  trees=100, #caret does built in 500 trees. ntrees doesn't matter.
  min_n=tune()
) %>%
  set_engine("ranger")

rf_grid<-grid_regular(mtry(range=c(1,30)),min_n(range=c(1,200)),levels=2)

#SVM ----

#svm with linear kernel
genome_svm<-svm_poly(
  mode="classification",
  cost=tune(),
  degree=1
) %>%
  set_engine("kernlab")

svm_grid<-grid_regular(cost(range=c(5,10)),
                       levels=2,
                       original = T)


#workflow to assess all models

model_list <- list(genome_lr, genome_rf, genome_svm)
hyperparam_list <- list(lr_grid, rf_grid, svm_grid)

tuned_models <- map2(.x = model_list,
                     .y = hyperparam_list, 
                     .f = model_tune,
                     recipe = std_recipe,
                     train_data = genome_train,
                     cv_fold = 5)

saveRDS(tuned_models, file = "./results/models_tuned.rds")
