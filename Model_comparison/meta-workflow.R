#This is a workflow to compare different ML models in different demographic scenarios

#Load libraries. pacman doesn't like tidymodels

pacman::p_load(tidyverse,vip)
library(tidymodels)
source("./Model_comparison/model_tune.R")
source("./Model_comparison/model_vip.R")

#load data from cleaning script

genomes = read_csv("./data/bt_cpop.csv")

#truncate dataset to contain response and predictors only. Used for model fitting.

genomes_SS = genomes %>% 
  dplyr::select(sweep,H_1:h123_11)
skimr::skim(genomes_SS)

#Partition dataset into training and testing sets.
set.seed(1066)
genome_split<-initial_split(genomes_SS,prop=0.8)
genome_train = training (genome_split)
genome_test = testing (genome_split)

#Recipe design ----

#The standard recipe just standardizes all the predictors (mean=0, var =1)

std_recipe <- recipe(sweep ~., data=genome_train) %>% #set sweep as response variable. everything else is a predictor.
  step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
  step_normalize(all_predictors()) %>% #normalize all predictors
  prep()

baked_train <- bake(std_recipe, genome_train)

#Designate model ----

#Logistical regression with L2 regularisation
genome_lr = logistic_reg(
  mode="classification",
  penalty = tune(),
  mixture=1
) %>%
  set_engine("glmnet")

#Hyperparameter tuning----

#Create set of tuning parameters
lr_grid = grid_regular(penalty(range=c(0,0.2)),
                       levels=100, 
                       original = F)

lr_results = model_tune(recipe = std_recipe, 
                        train_data = genome_train, 
                        cv_folds = 10, 
                        model = genome_lr , 
                        tuning_params = lr_grid, 
                        seed = 1)

lr_imp = model_vip(model = genome_lr, features = colnames(baked_train))




