#assess how well ML predicts bottleneck data when only trained on constant pop. simulations

pacman::p_load(tidyverse)
library(tidymodels)
library(vip)

genomes = read_csv("~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")

#data cleaning
genomes$sweep = as.factor(genomes$sweep)
levels(genomes$sweep) = c(1,0)

dt = data.table::data.table(genomes)
dt[, category :=.GRP, by = c("bottle_time1", "bottle_size1",
                             "bottle_time2", "bottle_size2")][]

dt[, category:= as.factor(category)]

#take our constant population simulations for training
cpop_data = subset(dt, category %in% "1")
#remove variables that shouldn't be predictors
cpop_data [, c("X1","ID","category"):=NULL] 
cpop_data[,(2:6):= NULL]
cpop_data = as.data.frame(cpop_data)

set.seed(1)
cpop_split = initial_split(cpop_data, prop = 0.8)
cpop_train = training(cpop_split)
cpop_test = testing (cpop_split)

#preprocessing

genomes_recipe <- recipe(sweep ~., data=cpop_train) %>%
  step_corr(all_predictors(),threshold = 0.8) %>%
  step_center(starts_with("H"), starts_with("D"), starts_with("L"), 
              starts_with("w"), starts_with("Z")) %>%
  step_scale(starts_with("H"), starts_with("D"), starts_with("L"), 
             starts_with("w"), starts_with("Z")) %>%
  prep()

#Logistical Regression with L2 regularisation ----

genome_lr = logistic_reg(
  mode="classification",
  penalty = tune(),
  mixture=1
) %>%
  set_engine("glmnet")

lr_grid = grid_regular(penalty(range=c(0,0.2)),
                       levels=100, 
                       original = F)

wkfl = workflow() %>%
  add_recipe(genomes_recipe)

wkfl_lr= wkfl %>% add_model(genome_lr)
cv_splits<-vfold_cv(cpop_train,v=10,strata="sweep")

lr_res = tune_grid(wkfl_lr,
                   resamples = cv_splits,
                   grid = lr_grid, 
                   metrics=metric_set(accuracy),
                   control=control_grid(save_pred = TRUE))

#cv accuracy of each LR model
lr_metrics = collect_metrics(lr_res)

#plot cv accuracy as a function of tuning params
ggplot(data= lr_metrics,
       aes(x=penalty, y=mean)) +
  geom_point()

#select tuning param with highest cv accuracy
best_tuning = lr_res %>% 
  select_best(metric="accuracy")

wkfl_best = finalize_workflow(wkfl_lr,best_tuning)

final_lr = fit(wkfl_best, data = cpop_train )

# Testing model on all the different population bottlenecks

#conf_mat(data= cpop_test, final_lr)

lr_pred = predict(final_lr, cpop_test)
truth<-cpop_test$sweep
confusionMatrix(lr_pred$.pred_class,cpop_test$sweep)

# temp = cbind(truth,lr_pred)
# roc_auc(.pred_1,data=temp)
# 
# temp %>% 
# roc_auc(truth=truth,estimate = lr_pred$.pred_1)
