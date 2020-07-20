#This is a workflow to compare different ML models in different demographic scenarios

#Load libraries. pacman doesn't like tidymodels

library(tidyverse)
library(vip)
library(tidymodels)
source("./Model_comparison/model_tune.R")
source("./Model_comparison/model_vip.R")
source("./Model_comparison/model_performance.R")

#load data from cleaning script

genomes = read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/split_snp/0.25ds_snp_set.csv")

#truncate dataset to contain response and predictors only. Used for model fitting.

genomes_SS = genomes %>% 
  dplyr::select(sweep,H_1:h123_11,demography,severity)
skimr::skim(genomes_SS)
genomes_SS$demography <- as.factor(genomes_SS$demography)
#genomes_SS$severity <- as.factor(genomes_SS$severity)


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
  update_role(severity, new_role = 'demography') %>% #remove demography as a predictor
  step_corr(all_predictors(),threshold = 0.9) %>% #remove all highly correlated predictors
  step_normalize(all_predictors()) %>% #normalize all predictors
  prep()

baked_train <- bake(std_recipe, genome_train)

#downsample the training dataset for VIP calculations
set.seed(2261941)
k = dim(baked_train)[1]*0.1
ds_baked_train = sample_n(baked_train, k)


#Designate model and hyperparameters ----

#Logistical regression with L1 regularisation ----
genome_lr = logistic_reg(
  mode="classification",
  penalty = tune(),
  mixture= 1
) %>%
  set_engine("glmnet")

#Create set of tuning parameters
lr_grid = grid_regular(penalty(range=c(0,0.1)) ,
                       levels=10, 
                       original = F)

lr_results = model_tune(recipe = std_recipe,
                        train_data = genome_train,
                        cv_folds = 10,
                        model = genome_lr ,
                        tuning_params = lr_grid,
                        seed = 162)

lr_imp = model_vip(model = lr_results, baked_data = ds_baked_train)
saveRDS(lr_imp, file = './results/lr_imp.rds')

#RandomForest classifier

genome_rf<-rand_forest(
  mode="classification",
  mtry=tune(),
  trees=500, #caret does built in 500 trees. ntrees doesn't matter.
  min_n=tune()
) %>%
  set_engine("ranger")

rf_grid<-grid_regular(mtry(range=c(10,80)),min_n(range=c(100,1000)),levels=4)

rf_results = model_tune(recipe = std_recipe,
                        train_data = genome_train,
                        cv_folds = 10,
                        model = genome_rf ,
                        tuning_params = rf_grid,
                        seed = 2)

rf_imp = model_vip(model = rf_results, baked_data = baked_train)
saveRDS(rf_imp, file = './results/rf_imp.rds')

#SVM ----

#svm with linear kernel
# genome_svm<-svm_poly(
#   mode="classification",
#   cost=tune(),
#   degree=1
# ) %>%
#   set_engine("kernlab")
# 
# svm_grid<-grid_regular(cost(range=c(10,20)),
#                        levels=1,
#                        original = T)
# 
# 
# svm_results = model_tune(recipe = std_recipe,
#                         train_data = genome_train,
#                         cv_folds = 10,
#                         model = genome_svm ,
#                         tuning_params = svm_grid,
#                         seed = 1)

#MARS ----
genome_mars <- mars(
  mode = "classification",
  prod_degree = tune(),
  num_terms = tune(),
  prune_method = "backward" #find default
) %>% 
  set_engine("earth")

n = 5
mars_grid = grid_regular(num_terms(range=c(1,110)), levels = n) %>%
  cbind(prod_degree = c(rep(1,n),rep(2,n)))

mars_results = model_tune(recipe = std_recipe,
                         train_data = genome_train,
                         cv_folds = 10,
                         model = genome_mars ,
                         tuning_params = mars_grid,
                         seed = 1)

mars_imp = model_vip(model = mars_results, baked_data = baked_train)

#RDA ----
genome_rda <- discrim::discrim_regularized(
  mode = 'classification', 
  frac_common_cov = tune(), #lambda
  frac_identity = tune() #gamma
) %>%
  set_engine("klaR")

#ref https://rdrr.io/cran/klaR/man/rda.html, https://discrim.tidymodels.org/reference/discrim_regularized.html

rda_grid <- grid_regular(frac_common_cov=discrim::frac_common_cov(range=c(0,1)),
                   discrim::frac_identity(range=c(0,1)),
                   levels=10)

#hack for weird bug
names(rda_grid)[1] <- "frac_identity"

rda_results <- model_tune(recipe = std_recipe,
                          train_data = genome_train,
                          cv_folds = 10,
                          model = genome_rda ,
                          tuning_params = rda_grid,
                          seed = 1)

#Workflow to assess all models----

model_list <- list(genome_lr, genome_rf, genome_mars, genome_rda)
hyperparam_list <- list(lr_grid, rf_grid, mars_grid, rda_grid)

doParallel::registerDoParallel()
#tune all the models
tuned_models <- map2(.x = model_list,
             .y = hyperparam_list, 
             .f = model_tune,
             recipe = std_recipe,
             train_data = genome_train,
             cv_fold = 10)


saveRDS(tuned_models, file = "./results/models_tuned.rds")

#find variables of importance for each model
vip_all <- map(.x = tuned_models,
               .f = model_vip,
               baked_data=ds_baked_train)

saveRDS(vip_all, file = "./results/vip_all.rds")
#11:27

#extract the finalised workflows for each model
finalised_models <- list()
for(i in 1:length(model_list)){
  finalised_models[[i]] <- tuned_models[[i]]$fitted_model
}

#check the performance of each model by mapping the model_performance()
model_robustness <- map(.x = finalised_models, 
                        .f = model_performance,
                        test_data = genome_test, 
                        recipe = std_recipe)

#attach names for each AUC tibble
for( i in 1:length(model_list)){
  #add names to each list
  names(model_robustness)[i] <- tuned_models[[i]]$model
  
  #add the ML method used for each AUC tibble
  model_robustness[[i]] <- model_robustness[[i]] %>%
    mutate(method = names(model_robustness)[i])
}

#bind all the AUC tibbles into the one dataframe
robustness_df <- do.call(rbind, model_robustness)

#auc plot across bottleneck severities
ggplot(data = robustness_df,
       aes(x = severity+1, y = .estimate, color = method)) + #+1 to offset severity 0
  geom_line() +
  scale_x_log10() + 
  ylab("AUC") +
  xlab("severity")

#plot of tuning params vs cv accuracy
lr_tune <- tuned_models[[1]]$tune_tibble

ggplot(data = lr_tune, 
       aes(x = penalty, y = mean, color = mixture)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err))

rf_tune <- tuned_models[[2]]$tune_tibble

ggplot(data = rf_tune, 
       aes(x = mtry, y = mean, color = min_n)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err)) +
  ylab("cv accuracy")

mars_tune <- tuned_models[[3]]$tune_tibble

ggplot( data = mars_tune,
        aes( x = num_terms, y = mean, color = prod_degree)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err))

rda_tune <- tuned_models[[4]]$tune_tibble

ggplot( data = rda_tune,
        aes( x = frac_common_cov, y = mean, color = frac_identity)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err))

#VIP for each model
vip_df <- list()
for (i in 1:length(vip_all)){
  vip_df[i] <- vip_all[[i]][2]
  
  names(vip_df)[i] <- tuned_models[[i]]$model
  
  vip_df[[i]] <- vip_df[[i]] %>%
    mutate(model = tuned_models[[i]]$model) %>% 
    mutate(SS_type = case_when(stringr::str_extract(Variable, "^.{2}") %in% c("H_", "D_") ~ 'SFS',
                               stringr::str_extract(Variable, "^.{2}") %in% c("h1", "h2") ~ 'Haplotype',
                               stringr::str_extract(Variable, "^.{3}") %in% c("Zns", "LD_","w_m") ~ 'LD',
                               TRUE ~ 'Other'))
}

vip_df <- do.call(rbind, vip_df)


# temp <- vip_df[[1]]
# temp %>% 
#   filter(str_detect(Variable, "D_"))
# 
# str_detect(name, "TajD")

vip_df$Variable <- vip_df$Variable %>% as.factor()

#plot SFS importance

vip_df %>% 
  dplyr::filter(SS_type == 'SFS') %>%
  ggplot ( aes(x = Variable, y = Importance, color = model)) +
  geom_point()

vip_df %>% 
  dplyr::filter(SS_type == 'Haplotype') %>%
  ggplot ( aes(x = Variable, y = Importance, color = model)) +
  geom_point()

vip_df %>% 
  dplyr::filter(SS_type == 'LD') %>%
  ggplot ( aes(x = Variable, y = Importance, color = model)) +
  geom_point()

##Hub's version

library(data.table)
vip_df <-data.table(vip_df)
gsub(x = vip_df$Variable, pattern = "_\\d+$", replacement = "")
vip_df[, Stat := gsub(x = Variable, pattern = "_\\d+$", replacement = "")]
vip_df[, Window := as.numeric(gsub(x = Variable, pattern = ".*_(\\d+$)", replacement = "\\1"))]

vip_df %>%
  dplyr::filter(SS_type == 'SFS') %>%
  ggplot(aes(x = Window, y = Importance, col = model)) +
  geom_point() +
  facet_wrap("Stat")

vip_df %>%
  dplyr::filter(SS_type == 'Haplotype') %>%
  ggplot(aes(x = Window, y = Importance, col = model)) +
  geom_point() +
  facet_wrap("Stat")

vip_df %>%
  dplyr::filter(Stat == 'w_max') %>%
  ggplot(aes(x = Window, y = Importance, col = model)) +
  geom_point() +
  facet_wrap("Stat")


#contrived shenanagens




# for( i in 1:length(model_list)){
#   #add names to each list
#   names(model_robustness)[i] <- tuned_models[[i]]$model
#   
#   #add the ML method used for each AUC tibble
#   model_robustness[[i]] <- model_robustness[[i]] %>%
#     mutate(method = names(model_robustness)[i])
# }

#grid.arrange(grobs = pdps, ncol = 4)

