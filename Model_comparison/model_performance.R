#' model_performance function
#' 
#' Takes a fitted model and tests it on data from a bunch of different demographies.
#' Outputs the AUC for each demographic scenario. 
#'
#' @param fitted_model : A parsnip fitted model for detecting sweeps
#' @param test_data : A dataframe of testing data with different demographies
#' @param recipe: Recipe used in the model fitting workflow.
#'
#' @return A tibble for the AUC under each demographic scenario
#' @export
#'
#' @examples
model_performance <- function (fitted_model, test_data, recipe){
  
  demographies <- unique(test_data$severity) %>% 
    sort()
  dem_auc=list()
  
  i = 1
  for (fac in demographies){
    dem_data <- test_data %>%
      dplyr::filter(severity == fac)
    preds <- predict(fitted_model, dem_data, type = 'prob')
    truth <- dem_data$sweep %>% as.factor()
    dem_auc[[i]] <-roc_auc(tibble(preds,truth), truth = truth, .pred_hard) %>%
      mutate(severity = fac)
    i = i+1
  }
  output_auc = do.call(rbind,dem_auc)
  return(output_auc)
}