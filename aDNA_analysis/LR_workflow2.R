#This is a workflow to assess how well LR can detect sweeps for different 
#aDNA damage params

library(tidyverse)
library(vip)
library(tidymodels)

source("./Model_comparison/model_tune.R")

#read in data
#ancient_genomes <- read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)

ancient_genomes$impute_method <- as.factor(ancient_genomes$impute_method)
random_imp_genomes = ancient_genomes %>%
  filter(impute_method=="random")

#skimr::skim(ancient_genomes)


#Partition dataset into training and testing sets.
set.seed(820)
genome_split<-initial_split(random_imp_genomes,prop=0.8)
randimp_genome_train = training (genome_split)
randimp_genome_test = testing (genome_split)

#Recipe

#Organise variables into groups for the recipe

#Since the haplotype statistics are proportions, we don't want to normalise them.
hap_cols <- colnames(ancient_genomes)[17:35]
aDNA_dmg_cols <- colnames(ancient_genomes)[4:6]

randimp_std_recipe <- recipe(sweep ~., data = randimp_genome_train) %>% #set sweep as response variable. everything else is a predictor.
  update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
  update_role(impute_method, new_role = 'damage') %>%
  update_role(ID, new_role = "ID") %>%
  update_role( aDNA_dmg_cols, new_role = 'damage') %>% 
  add_role(hap_cols, new_role = 'haplotype') %>%
  step_corr(all_predictors(),threshold = 0.9) %>% #remove all highly correlated predictors
  step_normalize(all_predictors(), -has_role("haplotype")) %>% #normalize all predictors, except haplotype stats
  prep()

#if you get Error in lognet, check if you have NAs in final transformed data

check <- summary(randimp_std_recipe)

#transform training data for vip functions
randimp_baked_train <- bake(randimp_std_recipe, randimp_genome_train)

#Model fitting
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

#fit model
randimp_lr_results = model_tune(recipe = randimp_std_recipe,
                        train_data = randimp_genome_train,
                        cv_folds = 5,
                        model = genome_lr ,
                        tuning_params = lr_grid,
                        seed = 162)

#AUC for each combination of s_coef and missing rate
source("./Model_comparison/auc_scoef.R")
auc_scoef(fitted_model = randimp_lr_results$fitted_model, 
          test_data = randimp_genome_test,
          recipe = randimp_std_recipe)

source("./aDNA_analysis/auc_aDNA.R")

randimp_auc_res = auc_aDNA(fitted_model = randimp_lr_results$fitted_model, 
                   test_data = randimp_genome_test,
                   recipe = randimp_std_recipe)

ggplot(data = randimp_auc_res, aes(x = missing_rate, y = .estimate, color = s_coef)) +
  geom_point() +
  ylab("AUC") +
  ggtitle("Random Imputation Deam")

source("./Model_comparison/model_vip.R")
lr_imp = model_vip(model = lr_results, baked_data = baked_train)

saveRDS(lr_imp, file = './results/lr_imp.rds')