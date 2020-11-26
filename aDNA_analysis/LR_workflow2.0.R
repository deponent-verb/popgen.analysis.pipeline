pacman::p_load(tidyverse,vip,drake,glmnet)
library(tidymodels)

#read in data
#ancient_genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes = read_csv("./data/cleaned_aDNA.csv")
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)
#ancient_genomes$ID <- as.factor(ancient_genomes$ID)
ancient_genomes = ancient_genomes %>%
  select(-ID)

#read in functions. Must not be done outside because drake calls functions from
#external environment. 
source("./aDNA_analysis/auc_aDNA.R")


genomes_group <- function(){
  #go make a new factor of combined factor variables
  ancient_genomes %>%
    mutate(tech = interaction(impute_method,denoise_method)) %>%
    split(f = .$tech)
}

fit_model <- function(genomes_group){
  genomes = genomes_group[[1]]
  set.seed(12)
  #split data into training and test set
  genome_split <- initial_split(genomes,prop=0.8,strata = missing_rate)
  genome_train = training(genome_split)
  genome_test = testing(genome_split)
  
  cv_splits = vfold_cv(genome_train, v = 4)
  
  #Since the haplotype statistics are proportions, we don't want to normalise them.
  hap_cols <- colnames(genomes)[which(colnames(genomes)=="h1_1"):which(colnames(genomes)=="h123_5")]
  aDNA_dmg_cols <- c(colnames(genomes)[which(colnames(genomes)=="missing_rate"):which(colnames(genomes)=="denoise_method")],"tech")
  
  
  std_recipe <- recipe(sweep ~., data = genome_train) %>% #set sweep as response variable. everything else is a predictor.
    update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
    #update_role(ID, new_role = "sim_ID") %>%
    update_role( all_of(aDNA_dmg_cols), new_role = 'damage') %>% 
    add_role(all_of(hap_cols), new_role = 'haplotype') %>%
    step_corr(all_predictors(),-has_role("sim_ID"),threshold = 0.9) %>% #remove all highly correlated predictors
    step_normalize(all_predictors(), -has_role("haplotype")) %>% #normalize all predictors, except haplotype stats
    prep()
  
  #Model fitting
  model = logistic_reg(
    mode="classification",
    penalty = tune(),
    mixture= 1
  ) %>%
    set_engine("glmnet")
  
  #Create set of tuning parameters
  tuning_grid = grid_regular(penalty(range=c(0,0.1)) ,
                         levels=2,
                         original = F)
  
  #setup workflow
  meta_workflow <- workflows::workflow() %>%
    workflows::add_recipe(std_recipe) %>%
    workflows::add_model(model)
  
  tuning = tune::tune_grid(meta_workflow,
                           resamples = cv_splits,
                           grid = tuning_grid,
                           metrics= yardstick::metric_set(accuracy),
                           #save out of sample predictions 
                           control=tune::control_grid(save_pred = TRUE))
  
  #get tibble of CV accuracy for the different tuning params
  tune_results = tune::collect_metrics(tuning)
  
  #find best tuning parameters
  best_params <- tuning %>%
    tune::select_best(metric = "accuracy")
  
  #fit final model on training data using the best set of tuning params
  final_workflow <- tune::finalize_workflow(meta_workflow, best_params) %>%
    parsnip::fit(data = genome_train)
  
  #compute AUC for each set of missing rates and s_coef
  tech = genomes$tech %>% unique()
  auc_res = auc_aDNA(fitted_model = final_workflow, 
                          test_data = genome_test,
                          recipe = std_recipe) %>%
    mutate(tech = tech)
}

plan <- drake_plan(
  genomes = genomes_group(),
  model = target(fit_model(genomes), dynamic = map(genomes))
)

make(plan)
readd(model)

results = readd(model)

#generate AUC tables

table = results %>%
  split(f =.$tech)

sapply(table, function(x){mean(x$.estimate)})

#generate AUC plots

ggplot(data = results, aes(x = missing_rate, y = .estimate, color = factor(s_coef))) +
  geom_point() +
  geom_line(aes(group = factor(s_coef))) +
  ylab("AUC") +
  facet_wrap("tech", scales = "free_y") +
  scale_color_discrete(name = "technique")

results_zero = results %>% 
  filter( (tech == "zero.cluster") | (tech == "zero.fixed_cluster") | (tech =="zero.none"))

ggplot(data = results_zero, aes(x = missing_rate, y = .estimate, color = tech)) +
  geom_point() +
  geom_line(aes(group = tech)) +
  ylab("AUC") +
  facet_wrap(~s_coef)+
  scale_color_discrete(name = "technique")

results_random = results %>% 
  filter( (tech == "random.cluster") | (tech == "random.fixed_cluster") | (tech =="random.none"))

ggplot(data = results_random, aes(x = missing_rate, y = .estimate, color = tech)) +
  geom_point() +
  geom_line(aes(group = tech)) +
  ylab("AUC") +
  facet_wrap(~s_coef)+
  scale_color_discrete(name = "technique")


