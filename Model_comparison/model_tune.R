#' model tune function
#' 
#' Fits a parsnip model on some input data, using a set of tuning parameters. Returns
#' the best model by cv accuracy and a tibble of cv accuracy for each set of tuning 
#' params. (update to add computational time)
#'
#' @param recipe : A recipe object from recipes. 
#' @param train_data : A dataframe of training data, with response variable sweep and pop gen summary statistics as predictors. 
#' @param model : A parsnip model object. 
#' @param tuning_params : A set of tuning parameters for the designated model.
#' @param cv_folds: Number of folds to break training data for tuning.
#' @param seed: optional random seed for cross validation
#'
#' @return List containing the model type (e.g. randomForest), the best set of tuning params,
#' a tibble of tuning results, computational time taken for tuning, a final worfklow with the 
#' model fitted on the whole training data using the best set of tuning params. 
#' @export
#'
#' @examples
model_tune <- function (recipe, train_data , model , tuning_params, cv_folds, seed = NA){

  #check inputs (not complete)  
  if(is.na(seed)) {
    seed = runif(1,0, 9999)
  }
  
  #Break training data into folds
  set.seed(seed)
  cv_splits<-vfold_cv(train_data,v=10,strata="sweep")
  
  # model workflow
  meta_workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(model)
  
  #unparallel version. prelim profiling shows parallelisation gives no speed bonus. Makes collecting metrics later more fiddly as well.
  
  t1 = Sys.time()
  tuning = tune_grid(meta_workflow,
                     resamples = cv_splits,
                     grid = tuning_params,
                     metrics=metric_set(accuracy),
                     control=control_grid(save_pred = TRUE))
  t2 = Sys.time()

  tune_results = collect_metrics(tuning)

  
  #find the best tuning parameters based on cv accuracy
  best_params <- tuning %>%
    select_best(metric = "accuracy")
  
  #finalize the workflow. Take model with the best set of tuning params and fit onto the whole training data. 
  
  final_workflow <- finalize_workflow(meta_workflow, best_params) %>%
    fit(data = train_data)
  
  #gather outputs
  source("./Model_comparison/tuned_model.R")
  output <- tuned_model(model = class(model)[1], 
                 best_params = best_params,
                 tune_tibble = tune_results, 
                 comp_time = t2-t1,
                 fitted_model = final_workflow)
  return(output)
}

# #set up parallel computing for tuning
# cores = parallel::detectCores()
# cl = parallel::makeCluster(cores)
# `%dopar%` <- foreach::`%dopar%`
# 
# #hyperparam tuning via grid search
# doParallel::registerDoParallel(cl,cores = cores)
# tune_blocks = split(tuning_params, as.factor(1:cores))
# 
# tuning_list = foreach::foreach (i = 1:length(tune_blocks), .packages = c('tune','yardstick')) %dopar% {
#   
#   tune::tune_grid(meta_workflow,
#             resamples = cv_splits,
#             grid = tune_blocks[[i]], 
#             metrics=metric_set(accuracy),
#             control=control_grid(save_pred = TRUE)) 
# }
# 
# #create a tibble for the cv accuracy of each set of tuning params
# tuning = lapply(tuning_list,collect_metrics) 
# tune_results = do.call(rbind, tuning)
# 

