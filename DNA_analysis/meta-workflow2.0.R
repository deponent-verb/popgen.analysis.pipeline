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
# genome_test = testing (genome_split)

cv_splits<-rsample::vfold_cv(genome_train,v=10,strata="sweep")

#The standard recipe just standardizes all the predictors (mean=0, var =1) 

std_recipe <- recipe(sweep ~., data = genome_train) %>% #set sweep as response variable. everything else is a predictor.
  update_role(demography, new_role = 'demography') %>% #remove demography as a predictor
  update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
  update_role(severity, new_role = 'demography') %>% #remove severity as a predictor
  step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
  step_normalize(all_predictors()) %>% #normalize all predictors
  prep()

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

#3.329075 mins for 10 levels

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

features = c("w_max_5","w_max_3","w_max_7","w_max_8")
lr_pdp = list()

for(f in 1:length(features)){
  lr_pdp[[f]] = pdp::partial(lr_caret, pred.var = features[f],
                             plot = TRUE, type = "classification")
}

grid.arrange(grobs = lr_pdp)

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
                         levels=10)
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

#1.04 hours for 100
