#This is a workflow to compare different ML models in different demographic scenarios

#Load libraries. pacman doesn't like tidymodels

pacman::p_load(tidyverse,vip)
library(tidymodels)
source("./Model_comparison/model_tune.R")
source("./Model_comparison/model_vip.R")
source("./Model_comparison/model_performance.R")

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

#remove demography as a predictor
# genome_train <- genome_train %>% 
#   select(-demography)
#update_role(sample, new_role = "id variable") %>%

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

# rf_results = model_tune(recipe = std_recipe, 
#                         train_data = genome_train, 
#                         cv_folds = 10, 
#                         model = genome_rf , 
#                         tuning_params = rf_grid, 
#                         seed = 1)
# 
# rf_imp = model_vip(model = rf_results, baked_data = baked_train)

#SVM ----

#svm with linear kernel
genome_svm<-svm_poly(
  mode="classification",
  cost=tune(),
  degree=1
) %>%
  set_engine("kernlab")

svm_grid<-grid_regular(cost(range=c(5,10)),
                       levels=6,
                       original = T)


# svm_results = model_tune(recipe = std_recipe, 
#                         train_data = genome_train, 
#                         cv_folds = 10, 
#                         model = genome_svm , 
#                         tuning_params = svm_grid, 
#                         seed = 1)

#MARS
# genome_mars <- mars(
#   mode = "classification",
#   prod_degree = 1, 
#   prune_method = default #find default
# )



#Running workflow functions on each model

#logistic regression example
# lr_results = model_tune(recipe = std_recipe, 
#                         train_data = genome_train, 
#                         cv_folds = 10, 
#                         model = genome_lr , 
#                         tuning_params = lr_grid, 
#                         seed = 1)
# 
# lr_auc = model_performance(fitted_model = lr_results$fitted_model, 
#                            test_data = genome_test,
#                            recipe = std_recipe)

# lr_imp = model_vip(model = lr_results, baked_data = baked_train)

#workflow to assess all models

model_list <- list(genome_lr, genome_rf, genome_svm)
hyperparam_list <- list(lr_grid, rf_grid, svm_grid)

tuned_models <- map2(.x = model_list,
             .y = hyperparam_list, 
             .f = model_tune,
             recipe = std_recipe,
             train_data = genome_train,
             cv_fold = 10)

saveRDS(tuned_models, file = "./results/models_tuned.rds")

#model_performance <- function (fitted_model, test_data, recipe)

##ignore below




