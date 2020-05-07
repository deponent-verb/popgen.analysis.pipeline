#' model_vip function
#' 
#' Takes a fitted model and the preprocessed data used for fitting. Plots the partial dependence plots
#' , ice plots and importance scores for each predictor in the fitted model.  
#'
#' @param model : A tuned_model object.
#' @param baked_data : A dataframe that has been transformed in a recipe. 
#'
#' @return Returns plots in a list
#' @export
#'
#' @examples
model_vip <- function(model, baked_data){
  
  if(class(model)!="tuned_model"){
    stop("argument model is not of class `tuned_model`")
  }
  
  #refit model in caret because vip functions are not compatible with parsnip models.
  model_type = model$model
  fitControl <- caret::trainControl(method = "none", classProbs = TRUE)
  
  if(model_type=="logistic_reg"){
    caret_model = caret::train(
      sweep~. ,
      data = baked_data,
      method = 'glmnet',
      trControl = fitControl,
      tuneGrid = data.frame(alpha = 1 ,lambda = model$best_params[[1]]),
      metric = "accuracy"
    )
  } else {
    stop("model type not supported.")
  }
  
  #partial dependence plots
  
  cols = colnames(baked_data)
  features = cols[!cols %in% "sweep"]
  pdps <- lapply(features, FUN = function(feature) {
    pd <- pdp::partial(caret_model, 
                       pred.var = feature, 
                       train = baked_train, 
                       type = "classification")
    autoplot(pd) +
      theme_light()
  })
  
  #ice curves
  
  ice_curves <- lapply(features, FUN = function(feature) {
    ice <- pdp::partial(caret_model, pred.var = feature, ice = TRUE)
    autoplot(ice, alpha = 0.1) + 
      theme_light()
  })
  
  #importance plot
  #imp_plot <- vip::vip(model, method = "firm", ice =T, train = baked_train)
  
  #importance measure
  imp <- vip::vi_firm(caret_model, feature_names = features)

  
  output <- list(pdps, ice_curves, imp)
  return(output)
}