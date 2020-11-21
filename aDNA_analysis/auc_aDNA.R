#' auc_aDNA function
#' 
#' Takes a fitted model and tests it on data from with different combinations of s_coef
#' and DNA damage parameters. Outputs the AUC for each selection coefficient and missing 
#' rate combination.
#'
#' @param fitted_model : A parsnip fitted model for detecting sweeps
#' @param test_data : A dataframe of testing data with different demographies
#' @param recipe: Recipe used in the model fitting workflow.
#'
#' @return A tibble for the AUC under each selection coefficient.
#' @export
#'
#' @examples
auc_aDNA <- function (fitted_model, test_data, recipe){
  
  #just in case tidyverse wasn't loaded
  `%>%` <- magrittr::`%>%`
  
  #take all selection coefficients
  s_values <- unique(test_data$s_coef) %>% 
    sort()
  missing_rates <- unique(test_data$missing_rate) %>%
    sort()
  #remove first element because it is 0. i.e. neutral. We test data on a combination of neutral and hard data.
  s_values <- s_values[-1]
  dem_auc=list()
  
  i = 1
  for (s in s_values){
    print(s)
    for (miss in missing_rates){
      neutral_data = test_data %>%
        dplyr::filter(s_coef == 0, missing_rate == miss)
      
      hard_data <- test_data %>%
        dplyr::filter(s_coef == s, missing_rate == miss)
      
      dem_data = dplyr::bind_rows(hard_data,neutral_data)
      truth <- dem_data$sweep %>% as.factor()
      preds <- predict(fitted_model, dem_data, type = 'prob')
      dem_auc[[i]] <-roc_auc(tibble(preds,truth), truth = truth, .pred_hard) %>%
        mutate(s_coef = s, missing_rate = miss)
      i = i+1
    }
  }
  output_auc = do.call(rbind,dem_auc)
  return(output_auc)
}