#Load libraries. pacman doesn't like tidymodels

pacman::p_load(tidyverse)
library(tidymodels)

#Read in data ----

genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/bt_cpop.csv")

#remove bottleneck and selection coefficient params. We shouldn't fit model on these variables. 

genomes_SS = genomes %>% select(sweep,H_1:h123_11)
genomes_SS$sweep = as.factor(genomes_SS$sweep)
skimr::skim(genomes_SS)

#Partition dataset
genome_split<-initial_split(genomes_SS,prop=0.8)

#Preprocessing ----

#Standardize H and D. We don't transform the haplotype stats as they are already between 0,1. 

genomes_recipe <- recipe(sweep ~., data=training(genome_split)) %>%
  step_corr(all_predictors(),threshold = 0.8) %>%
  step_center(starts_with("H"),starts_with("D")) %>%
  step_scale(starts_with("H"),starts_with("D")) %>%
  prep()

genomes_recipe

genome_training <- bake(genomes_recipe, training(genome_split))
genome_testing<-bake(genomes_recipe,testing(genome_split))

#Model Fitting ----

#Pick tuning params based on 10-fold CV

set.seed(1688)
cv_splits<-vfold_cv(training(genome_split),v=10,strata="sweep")

#Logistical Regression with L2 regularisation ----

genome_lr = logistic_reg(
  mode="classification",
  penalty = tune(),
  mixture=1
) %>%
  set_engine("glmnet")

lr_grid = grid_regular(penalty(range=c(0,50)),
                       levels=50, 
                       original = F)

lr_res = tune_grid(genomes_recipe,
                   model=genome_lr,
                   resamples = cv_splits,
                   grid=lr_grid)

collect_metrics(lr_res)

best_tuning = lr_res %>% 
  select_best(metric="accuracy")

final_lr = finalize_model(genome_lr, best_tuning) %>%
  fit(sweep~., data=genome_training)

# Random Forest----

rf_grid<-grid_regular(mtry(range=c(10,30)),min_n(range=c(20,40)),levels=5)

genome_rf<-rand_forest(
  mode="classification",
  mtry=tune(),
  trees=500,
  min_n=tune()
) %>%
  set_engine("randomForest")

doParallel::registerDoParallel()
a = Sys.time()
rf_res<-tune_grid(genomes_recipe,
                  model=genome_rf,
                  resamples=cv_splits,
                  grid=rf_grid)
b = Sys.time()

best_tuning = rf_res %>% 
  select_best(metric="accuracy")

final_rf = finalize_model(genome_rf, best_tuning) %>%
  fit(sweep~., data=genome_training)


#SVM ----

genome_svm<-svm_poly(
  mode="classification",
  cost=tune(),
  degree=tune()
) %>%
  set_engine("kernlab")

svm_grid<-grid_regular(cost(range=c(7,9)),
                       degree(),
                       levels=3,
                       original = T)

doParallel::registerDoParallel()
svm_res<-tune_grid(genomes_recipe,
                   model=genome_svm,
                   resamples=cv_splits,
                   grid=svm_grid)

best_tuning<-svm_res %>%
  select_best(metric="accuracy")

final_svm<-finalize_model(genome_svm,best_tuning) %>%
  fit(sweep~., data=genome_training)


# Model Assessment

#obtain soft classification on testing set
rf_pred<-predict(final_rf,genome_testing,type="prob")
lr_pred<-predict(final_lr,genome_testing,type="prob")
svm_pred<-predict(final_svm,genome_testing,type="prob")
model_preds = list(lr_pred, rf_pred, svm_pred )

model_names<-c("logistic regression","random forest","support vector machine")

pacman::p_load(cowplot)
#grab true responses in test set
truth<-genome_testing$sweep

model_eval<-function(model_name,model_pred,truth){ 
  data<-as.data.frame(model_pred)
  ans<-data %>% mutate(model=model_name) %>% cbind(truth) 
  return(ans)
}

model_list<-list()
for(i in 1:length(model_names)){
  model_list[[i]]<-model_eval(model_names[i],model_preds[i],truth) 
}
model_out<-plyr::ldply(model_list, data.frame)

model_out %>%
  #get individual roc curves for each model 
  group_by(model) %>% 
  roc_curve(truth=truth, estimate=.pred_hard) %>% 
  ggplot(
    aes(x=1-specificity,y=sensitivity,color=model) 
  )+
  geom_line(size=1.1) +
  geom_abline(slope=1, intercept = 0, size=0.2) +
  #fix aspect ratio to 1:1 for visualization coord_fixed() +
  theme_cowplot() +
  xlab("FPR") +
  ylab("TPR")

# AUC

model_auc <- model_out %>%
  group_by(model) %>%
  roc_auc(truth=truth,.pred_hard,options = list(smooth = TRUE))

knitr::kable(model_auc)
