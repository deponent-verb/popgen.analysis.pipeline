#' model_vip function
#' 
#' Takes a fitted model and the preprocessed data used for fitting. Plots the partial dependence plots
#' , ice plots and importance scores for each predictor in the fitted model.  
#'
#' @param model : A fitted model object from tidymodels
#' @param baked_data : A dataframe that has been transformed in a recipe. 
#'
#' @return Returns plots in a list
#' @export
#'
#' @examples
model_vip <- function(model, baked_data){
  
  #refit model in caret because vip functions are not compatible with parsnip models.
  model_fit <- model$fit
  model_type <- class(model_fit)
  
  #partial dependence plots
  
  pdps <- lapply(features, FUN = function(feature) {
    pd <- pdp::partial(model, 
                       pred.var = feature, 
                       train = baked_train, 
                       type = "classification")
    autoplot(pd) +
      theme_light()
  })
  
  #ice curves
  
  ice_curves <- lapply(features, FUN = function(feature) {
    ice <- pdp::partial(model, pred.var = feature, ice = TRUE)
    autoplot(ice, alpha = 0.1) + 
      theme_light()
  })
  
  #importance plot
  imp <- vip::vip(model, method = "firm", ice =T, train = baked_train)
  
  output <- list(pdps, ice_curves, imp)
  return(output)
}