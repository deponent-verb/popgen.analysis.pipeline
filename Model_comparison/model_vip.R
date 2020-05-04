#' model_vip function
#' 
#' Takes a model and a vector of features/predictors. Plots the partial dependence plots
#' , ice plots and importance scores for each feature. 
#'
#' @param model : A fitted model object from tidymodels
#' @param features : A vector of features for doing the plots
#'
#' @return Returns plots in a list
#' @export
#'
#' @examples
model_vip <- function(model, features){
  
  #partial dependence plots
  
  pdps <- lapply(features, FUN = function(feature) {
    pd <- pdp::partial(model, pred.var = feature)
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
  imp <- vip::vip(model, method = "firm", ice =T)
  
  output <- list(pdps, ice_curves, imp)
  return(output)
}