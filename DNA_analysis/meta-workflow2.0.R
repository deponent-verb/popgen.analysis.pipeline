#version 2 of workflow for chunkybit 1

pacman::p_load(tidyverse,vip)
library(tidymodels)

#load data from cleaning script

genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_chunky1_data.csv")
genomes = subset(genomes,
                 select = -c(ID))
genomes$sweep <- ifelse(genomes$sweep=="hard",1,0)
genomes$sweep <- as.factor(genomes$sweep)

#Partition dataset----
set.seed(1066)
genome_split<-initial_split(genomes,prop=0.8)
genome_train = training (genome_split)
genome_test = testing (genome_split)

cv_splits<-rsample::vfold_cv(genome_train,v=10,strata="sweep")

#we don't standardise haplotype stats because they are bounded by 0,1
hap_cols <- colnames(genomes)[which(colnames(genomes)=="h1_1"):which(colnames(genomes)=="h123_11")]

std_recipe <- recipe(sweep ~., data=genome_train) %>% #set sweep as response variable. everything else is a predictor.
  update_role(demography, new_role = 'demography') %>% #remove demography as a predictor
  update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
  update_role(severity, new_role = 'demography') %>% #remove severity as a predictor
  add_role(all_of(hap_cols), new_role = 'haplotype') %>%
  step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
  step_normalize(all_predictors(), -has_role("haplotype")) %>% #normalize all predictors, except haplotype stats
  prep()

#transform dataset for pdp

baked_data = bake(std_recipe, new_data = genome_train)
trans_data = baked_data[,which(colnames(baked_data)=="H_1"):which(colnames(baked_data)=="h2_11")] %>%
  cbind(sweep = baked_data$sweep)

#Logistical regression with L1 regularisation ----
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

doParallel::registerDoParallel()

lr_t1 = Sys.time()

lr_wkfl = workflow() %>%
  add_recipe(std_recipe) %>%
  add_model(genome_lr)

lr_tune = lr_wkfl %>%
  tune_grid(
    resamples = cv_splits,
    grid = lr_grid
  )
lr_t2 = Sys.time()

#3.860783 mins for 10 levels

#LR CV surface plots

lr_cv = lr_tune %>%
tune::collect_metrics() %>%
  dplyr::filter(.metric == "accuracy")


ggplot(data = lr_cv, 
       aes(x = penalty, y = mean)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err))+
  ylab("cv accuracy") +
  xlab("lambda") 

#finalise LR
best_lr <- lr_tune %>%
  tune::select_best(metric = "accuracy")

lr_final = lr_wkfl %>%
  tune::finalize_workflow(best_lr) %>%
  parsnip::fit(data = genome_train)

#vip

# lr_final %>%
#   pull_workflow_fit() %>%
#   vip(method = "firm", train = bake(std_recipe,genome_train), type = "classification")
#   
# lr_final %>%
#   pull_workflow_fit() %>%
#   vip(method = "firm", train = bake(std_recipe,genome_train),
#       type = "classification")
# 
# lr_fit = lr_final %>%
#   pull_workflow_fit()
# 
# vip(lr_fit, lamda = lr_fit$fit$lambda[100], method = "firm",
#     train = bake(std_recipe,genome_train), type = "classification",
#     new_data = NULL,
#     pred_wrapper = glmnet::predict.glmnet,
#     arg = "nonzero")

#not ideal. Inserted best param here. Couldn't get vip to work with glmnet.
caret_data = bake(std_recipe,genome_train)
caret_data = subset(caret_data,
                    select = -c(s_coef,demography,severity))
caret_data$sweep = ifelse(caret_data$sweep==1,"hard","neutral") %>% as.factor()

lr_caret = caret::train(
  sweep~. ,
  data = caret_data,
  method = 'glmnet',
  trControl = caret::trainControl(method = "none", classProbs = TRUE),
  tuneGrid = data.frame(alpha = 1 ,lambda = 0),
  metric = "accuracy"
)

#vip(lr_caret)

lr_imp = vip(lr_caret, method = "firm")

lr_imp +
  ggtitle("Logistic Regression")

#pdp

features = c("w_max_5","w_max_3","w_max_7")
lr_pdp = list()

for(f in 1:length(features)){
  lr_pdp[[f]] = pdp::partial(lr_caret, pred.var = features[f],
                             plot = TRUE, type = "classification")
}

grid.arrange(grobs = lr_pdp, ncol = 3)

#ice example
ice_t1 = Sys.time()
lr_ice = pdp::partial(lr_caret, pred.var = "w_max_5", 
                      ice = TRUE, type = "classification", plot = TRUE)
ice_t2 = Sys.time()
# ice_curves <- lapply(features, FUN = function(feature) {
#   ice <- pdp::partial(caret_model, pred.var = feature, ice = TRUE)
#   autoplot(ice, alpha = 0.1) + 
#     theme_light()
# })

# pdp::partial(lr_caret, pred.var = "w_max_5",
#              plot = TRUE, type = "classification")


#RDA ----
library("discrim")

genome_rda <- discrim::discrim_regularized(
  mode = 'classification', 
  frac_common_cov = tune(), #lambda
  frac_identity = tune() #gamma
) %>%
  set_engine("klaR")

#ref https://rdrr.io/cran/klaR/man/rda.html, https://discrim.tidymodels.org/reference/discrim_regularized.html

rda_grid <- grid_regular(frac_common_cov=discrim::frac_common_cov(range=c(0,1)),
                         discrim::frac_identity(range=c(0,1)),
                         levels=5)
names(rda_grid)[1] <- "frac_identity" #hack for weird bug

#fit rda
rda_wkfl = workflow() %>%
  add_recipe(std_recipe) %>%
  add_model(genome_rda)

rda_t1 = Sys.time()

rda_tune = rda_wkfl %>%
  tune_grid(
    resamples = cv_splits,
    grid = rda_grid
  )
rda_t2 = Sys.time()

#15.87438 mins for 25, 0.635 mins/model

rda_cv = rda_tune %>%
  tune::collect_metrics() %>%
  dplyr::filter(.metric == "accuracy")

ggplot(data = rda_cv,
       aes( x = frac_common_cov, y = mean, color = factor(frac_identity))) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err)) +
  xlab("lambda") +
  labs(color = "gamma") +
  ylab("cv accuracy")

#finalise rda
best_rda <- rda_tune %>%
  tune::select_best(metric = "accuracy")

rda_final = rda_wkfl %>%
  tune::finalize_workflow(best_rda) %>%
  parsnip::fit(data = genome_train)

#vip, caret version (not ideal)

trans_data1 = trans_data
trans_data1$sweep = ifelse(trans_data1$sweep==1,"hard","neutral")

rda_caret = caret::train(sweep~.,
  data = trans_data1,
  method = "rda",
  trControl = caret::trainControl(method = "none", classProbs = TRUE),
  tuneGrid = data.frame(gamma = best_rda$frac_identity,
                        lambda = best_rda$frac_common_cov))

rda_imp = vip(rda_caret, method = "firm")
saveRDS(rda_imp,file = "./results/Chunky1/rda_firm.rds")

rda_imp + 
  ggtitle("RDA")

features = c("D_6","H_6","D_7")
rda_pdp = list()

for(f in 1:length(features)){
  rda_pdp[[f]] = pdp::partial(rda_caret, pred.var = features[f],
                             plot = TRUE, type = "classification")
}

grid.arrange(grobs = rda_pdp, ncol = 3)

#vip

#does not work. Will have to use caret.

# rda_imp = rda_final %>%
#   pull_workflow_fit() %>%
#   vip(method = "firm", train = bake(std_recipe, genome_train),
#       features_names = c("w_max_6","w_max_5")) 
# 
# rda_final %>%
#   pull_workflow_fit() %>%
#   .$fit %>%
#   vip(method = "firm",train = bake(std_recipe, genome_train),
#       features_names = c("w_max_6","w_max_5"))

#Random Forest, test code with small hyperparams ---

#get number of predictors after applying recipe
num_terms = which(std_recipe$term_info$role == "predictor") %>% length()

genome_rf = rand_forest(
  mode = "classification",
  mtry = tune(),
  min_n = tune(),
  trees = 100 
) %>%
  set_engine("ranger")

rf_wkfl = workflow() %>%
  add_recipe(std_recipe) %>%
  add_model(genome_rf) 

rf_grid<-grid_regular(mtry(range=c(10,num_terms)),
                      min_n(range=c(100,1000)),levels=4)

rf_t1 = Sys.time()
rf_tune = rf_wkfl %>%
  tune_grid(
    resamples = cv_splits,
    grid = rf_grid
  )
rf_t2 = Sys.time()

#takes 0.2 hours for one model

rf_cv = rf_tune %>%
  collect_metrics() %>%
  dplyr::filter(.metric == "accuracy")

write_csv(rf_cv, file = "./results/Chunky1/rf_res.csv")

ggplot(data = rf_cv,
       aes( x = mtry, y = mean, color = cut(min_n, breaks = length(min_n) ))) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err)) +
  ylab("cv accuracy") + 
  labs(color = "min_n")

ggplot(data = rf_cv,
       aes( x = mtry, y = mean, color = factor(min_n))) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err)) +
  ylab("cv accuracy") + 
  labs(color = "min_n")

#finalise rf
best_rf <- rf_tune %>%
  tune::select_best(metric = "accuracy")

rf_final = rf_wkfl %>%
  tune::finalize_workflow(best_rf) %>%
  parsnip::fit(data = genome_train)

saveRDS(rf_final, file = "./results/Chunky1/rf_model.rds")

#vip
rf_imp = rf_final %>%
  pull_workflow_fit() %>%
  vip(method = "firm", target = "sweep", metric = "accuracy",
      pred_wrapper = ranger::predictions,
      train = bake(std_recipe,genome_train),
      new_data = NULL)

rf_imp +
  ggtitle("Random Forest")

saveRDS(rf_imp, file = "./results/Chunky1/rf_firm.rds")

features = c("D_6", "D_5", "D_4")
rf_pdp = list()

t1 = Sys.time()
for (f in 1:length(features)){
  rf_pdp[[f]] = rf_final %>%
    pull_workflow_fit() %>%
    .$fit %>%
    pdp::partial(train = bake(std_recipe,genome_train),
                 pred_wrapper = ranger::predictions,
                 pred.var = features[f],
                 plot = TRUE, 
                 type = "classification",
                 new_data = NULL,
                 which.class = 2)
}

grid.arrange(grobs = rf_pdp, ncol = 3)
t2 = Sys.time()

t1= Sys.time()
rf_final %>%
  pull_workflow_fit() %>%
  .$fit %>%
  pdp::partial(train = bake(std_recipe,genome_train),
               pred_wrapper = ranger::predictions,
               pred.var = "D_7",
               plot = TRUE, 
               type = "classification",
               new_data = NULL,
               which.class = 2)
t2 = Sys.time()


#MARS ----

genome_mars <- mars(
  mode = "classification",
  prod_degree = tune(),
  num_terms = tune(),
  prune_method = "forward" #find default
) %>% 
  set_engine("earth")

#get number of predictors after applying recipe
num_terms = which(std_recipe$term_info$role == "predictor") %>% length()

n = 5
mars_grid = grid_regular(num_terms(range=c(1,num_terms)), levels = n) %>%
  cbind(prod_degree = c(rep(1,n),rep(2,n)))

mars_t1 = Sys.time()

mars_wkfl = workflow() %>%
  add_recipe(std_recipe) %>%
  add_model(genome_mars)

mars_tune = mars_wkfl %>%
  tune_grid(
    resamples = cv_splits,
    grid = mars_grid
  )

mars_t2 = Sys.time()
#avg 0.847 mins for one model

#CV plots

mars_cv = mars_tune %>%
  collect_metrics() %>%
  dplyr::filter(.metric == "accuracy")

write_csv(mars_cv, file = "./results/Chunky1/mars_res.csv")

ggplot(data = mars_cv,
       aes( x = num_terms, y = mean, color = prod_degree)) +
  geom_point() +
  geom_errorbar( aes(ymax = mean + std_err, ymin = mean - std_err)) +
  ylab("cv accuracy") + 
  labs(color = "degree")

#finalise mars
best_mars <- mars_tune %>%
  tune::select_best(metric = "accuracy")

mars_final = mars_wkfl %>%
  tune::finalize_workflow(best_mars) %>%
  parsnip::fit(data = genome_train)

saveRDS(mars_final, file = "./results/Chunky1/mars_model.rds")

#mars vip

mars_imp = mars_final %>%
  pull_workflow_fit() %>%
  vip(method = "firm", train = bake(std_recipe,genome_train))

saveRDS(mars_imp, file = "./results/Chunky1/mars_firm.rds")
  
mars_imp + 
  ggtitle("MARS")

#mars pdp

features = c("H_6", "h1_6", "H_5")
mars_pdp = list()

for(f in 1:length(features)){
  mars_pdp[[f]] = mars_final %>%
    #pull parsnip model
    pull_workflow_fit() %>%
    #pull out MARS model since pdp does not have support for parsnip
    .$fit %>%
    pdp::partial(train = trans_data, pred.var = features[f], 
                 plot = TRUE, type = "classification")
}

grid.arrange(grobs = mars_pdp, ncol = 3)

ggplot(trans_data, aes(factor(sweep), D_6)) + 
  geom_boxplot()

#plot all firm scores
firm_plots = list(lr_imp,rda_imp,rf_imp,mars_imp)
grid.arrange(grobs = list(lr_imp + ggtitle("Logistic Regression"),
                          rda_imp + ggtitle("RDA"),
                          rf_imp + ggtitle("Random Forest"),
                          mars_imp + ggtitle("MARS")), 
             ncol = 2)


#overall AUC for each model

preds <- predict(list(lr_final,rda_final,rf_final,mars_final), genome_test, type = 'prob')
truth <- as.factor(genome_test$sweep)
roc_auc(tibble(preds,truth), truth = truth, .pred_0)

models = list(lr_final,rda_final,rf_final,mars_final)
  
for (m in models){
  preds <- predict(m, genome_test, type = 'prob')
  truth <- as.factor(genome_test$sweep)
  print(roc_auc(tibble(preds,truth), truth = truth, .pred_0)
)
}

#AUC for each severity
source("./Model_comparison/model_performance.R")

# model_performance(fitted_model = mars_final,
#                   test_data = genome_test,recipe = std_recipe)

models = list(lr_final,rda_final,rf_final,mars_final)
model_names = c("Logistic Regression", "RDA","Random Forest","MARS")
#check the performance of each model by mapping the model_performance()
model_robustness <- map(.x = models, 
                        .f = model_performance,
                        test_data = genome_test, 
                        recipe = std_recipe)

  
#attach names for each AUC tibble
for( i in 1:length(models)){
  #add names to each list
  names(model_robustness)[i] <- model_names[i]
  
  #add the ML method used for each AUC tibble
  model_robustness[[i]] <- model_robustness[[i]] %>%
    mutate(method = names(model_robustness)[i])
}

#bind all the AUC tibbles into the one dataframe
robustness_df <- do.call(rbind, model_robustness)

#auc plot across bottleneck severities
ggplot(data = robustness_df,
       aes(x = severity+1, y = .estimate, color = method)) + #+1 to offset severity 0
  geom_point() +
  scale_x_log10() + 
  ylab("AUC") +
  xlab("severity") + 
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))

#AUC for s_coef
source("./Model_comparison/auc_scoef.R")

model_scoef <- map(.x = models, 
                        .f = auc_scoef,
                        test_data = genome_test, 
                        recipe = std_recipe)


#attach names for each AUC tibble
for( i in 1:length(models)){
  #add names to each list
  names(model_scoef)[i] <- model_names[i]
  
  #add the ML method used for each AUC tibble
  model_scoef[[i]] <- model_scoef[[i]] %>%
    mutate(method = names(model_scoef)[i])
}

#bind all the AUC tibbles into the one dataframe
scoef_df <- do.call(rbind, model_scoef)

ggplot(data = scoef_df,
       aes(x = s_coef, y = .estimate, color = method)) +
  geom_point() +
  scale_x_log10() + 
  ylab("AUC") +
  xlab("selection coefficient") + 
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"))
