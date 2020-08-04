#' auc_scoef function
#' 
#' Takes a fitted model and tests it on data from a bunch of different selection coefficients.
#' The aim is to see how well the model distinguishes neutral data from sweeps of different 
#' strengths. Outputs the AUC for each selection coefficient. 
#'
#' @param fitted_model : A parsnip fitted model for detecting sweeps
#' @param test_data : A dataframe of testing data with different demographies
#' @param recipe: Recipe used in the model fitting workflow.
#'
#' @return A tibble for the AUC under each selection coefficient.
#' @export
#'
#' @examples
auc_scoef <- function (fitted_model, test_data, recipe){
  
  #take all selection coefficients
  s_values <- unique(test_data$s_coef) %>% 
    sort()
  #remove first element because it is 0. i.e. neutral
  s_values <- s_values[-1]
  dem_auc=list()
  
  neutral_data = test_data %>%
    dplyr::filter(s_coef == 0)
  
  i = 1
  for (s in s_values){
    print(s)
    hard_data <- test_data %>%
      dplyr::filter(s_coef == s)
    dem_data = dplyr::bind_rows(hard_data,neutral_data)
    truth <- dem_data$sweep %>% as.factor()
    preds <- predict(fitted_model, dem_data, type = 'prob')
    dem_auc[[i]] <-roc_auc(tibble(preds,truth), truth = truth, .pred_hard) %>%
      mutate(s_coef = s)
    i = i+1
  }
  output_auc = do.call(rbind,dem_auc)
  return(output_auc)
}