#script for MARS model fitting to ensure all analysis is done using the same
#fitting procedure

#Takes my ancient genome dataframe and fits a MARS model (output). 

#this is not a proper function!

MARS_fit <- function(genomes){
  set.seed(12)
  #split data into training and test set
  genome_split <- initial_split(genomes,prop=0.8,strata = missing_rate)
  genome_train = training(genome_split)
  genome_test = testing(genome_split)
  
  cv_splits = vfold_cv(genome_train, v = 10)
  
  #Since the haplotype statistics are proportions, we don't want to normalise them.
  hap_cols <- colnames(genomes)[which(colnames(genomes)=="h1_1"):which(colnames(genomes)=="h123_5")]
  aDNA_dmg_cols <- c(colnames(genomes)[which(colnames(genomes)=="missing_rate"):which(colnames(genomes)=="denoise_method")],"tech")
  
  
  std_recipe <- recipe(sweep ~., data = genome_train) %>% #set sweep as response variable. everything else is a predictor.
    update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
    #update_role(ID, new_role = "sim_ID") %>%
    update_role( all_of(aDNA_dmg_cols), new_role = 'damage') %>% 
    add_role(all_of(hap_cols), new_role = 'haplotype') %>%
    step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
    step_normalize(all_predictors(), -has_role("haplotype")) %>% #normalize all predictors, except haplotype stats
    prep()
  
  #Model fitting
  model <- mars(
    mode = "classification",
    prod_degree = 1,
    num_terms = tune(),
    prune_method = "backward" #find default
  ) %>% 
    set_engine("earth")
  
  #Create set of tuning parameters
  n = 5
  # tuning_grid = grid_regular(num_terms(range=c(2,10)), levels = n) %>%
  #   cbind(prod_degree = c(rep(1,n),rep(2,n)))
  tuning_grid = grid_regular(num_terms(range=c(4,15)), levels = n)
  
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
  
  return(final_workflow)
}

