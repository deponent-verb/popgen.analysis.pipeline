#Load libraries. pacman doesn't like tidymodels

pacman::p_load(tidyverse)
library(tidymodels)
library(vip)

#Read in data ----

genomes = read_csv("~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")

#remove bottleneck and selection coefficient params. We shouldn't fit model on these variables. 

genomes_SS = genomes %>% select(sweep,H_1:h123_11)
genomes_SS$sweep = as.factor(genomes_SS$sweep)
skimr::skim(genomes_SS)

#Partition dataset
set.seed(1707)
genome_split<-initial_split(genomes_SS,prop=0.8)
genome_train = training (genome_split)
genome_test = testing (genome_split)

#Preprocessing ----

#Standardize H and D. We don't transform the haplotype stats as they are already between 0,1. 

genomes_recipe <- recipe(sweep ~., data=genome_train) %>%
  step_corr(all_predictors(),threshold = 0.8) %>%
  step_center(starts_with("H"), starts_with("D"), starts_with("L"), 
              starts_with("w"), starts_with("Z")) %>%
  step_scale(starts_with("H"), starts_with("D"), starts_with("L"), 
             starts_with("w"), starts_with("Z")) %>%
  prep()

# genomes_recipe2 <- recipe (sweep ~., data= genome_train) %>%
#   step_corr(all_predictors(),threshold = 0.8) %>%
#   step_center(all_predictors()) %>%
#   step_scale(all_predictors()) %>%
#   prep()

genomes_recipe

# genome_training <- bake(genomes_recipe, training(genome_split))
# genome_testing<-bake(genomes_recipe,testing(genome_split))

#Model Fitting ----

#Pick tuning params based on 10-fold CV

set.seed(1688)
cv_splits<-vfold_cv(genome_train,v=10,strata="sweep")

#Logistical Regression with L2 regularisation ----

#workflows method

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

lr_start= Sys.time()
lr_res = tune_grid(wkfl_lr,
                   resamples = cv_splits,
                   grid = lr_grid, 
                   metrics=metric_set(accuracy),
                   control=control_grid(save_pred = TRUE))
lr_fin = Sys.time()

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

final_lr = fit(wkfl_best, data = genome_train )

#variables of importance

final_lr %>% 
  pull_workflow_fit() %>%
  tidy()

final_lr %>%
  pull_workflow_fit() %>%
  vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x=Importance, y = Variable, fill=Sign)) +
  geom_col() + 
  scale_x_continuous(expand = c(0,0)) + 
  labs(y=NULL) + 
  ggtitle("Logistic Regression VoI")
  


# Random Forest----

rf_grid<-grid_regular(mtry(range=c(1,30)),min_n(range=c(1,200)),levels=10)


genome_rf<-rand_forest(
  mode="classification",
  mtry=tune(),
  trees=50,
  min_n=tune()
) %>%
  set_engine("randomForest")

wkfl_rf = wkfl %>% add_model(genome_rf)

doParallel::registerDoParallel()
rf_start = Sys.time()

rf_res = tune_grid(wkfl_rf,
                  resamples=cv_splits,
                  grid=rf_grid,
                  metrics = metric_set(accuracy),
                  control=control_grid(save_pred = TRUE))

rf_end = Sys.time()

rf_metrics = collect_metrics(rf_res)

#plot cv accuracy as function of tuning params

rf_df = rf_metrics %>% 
  select(mtry,min_n,mean) 

rf_df$quart = cut(rf_df$mean, quantile(rf_df$mean), include.lowest = T)

ggplot(rf_df,aes(x = mtry, y = min_n, color = factor(quart))) +
  geom_point()

ggplot(rf_df,aes(x = mtry, y = min_n, z = mean)) +
  geom_density_2d()

#3d plot for better visualisation

x = rf_grid$mtry %>% unique()
y = rf_grid$min_n %>% unique()
n1 = length(x)
n2 = length(y)
cv = matrix(0, n1, n2)


for(i in 1:n1){
  for(j in 1:n2){
    print(i)
    print(j)
    row = rf_metrics %>% filter(mtry==x[i], min_n==y[j])
    cv[i,j] = row$mean
  }
}

fig = plot_ly(x = ~x, y = ~y,z = ~cv) %>% 
  add_surface()

fig = fig %>% layout(
  title = "CV accuracy as a function of RF tuning params",
  scene = list(
    xaxis = list(title = "mtry"),
    yaxis = list (title = "min_n")
  )
)

fig

#get best rf model
best_tuning = rf_res %>% 
  select_best(metric="accuracy")

wkfl_best_rf = finalize_workflow(wkfl_rf, best_tuning)

final_rf = fit( wkfl_best_rf, data= genome_train)

# RF VIP

final_rf %>% 
  pull_workflow_fit() %>%
  tidy()

final_rf %>%
  pull_workflow_fit() %>%
  vi() %>%
  mutate(
    Importance = abs(Importance),
    Variable = fct_reorder(Variable, Importance)
  ) %>%
  ggplot(aes(x=Importance, y = Variable)) +
  geom_col() + 
  scale_x_continuous(expand = c(0,0)) + 
  labs(y=NULL) + 
  ggtitle("RF VoI")


#SVM ----

genome_svm<-svm_poly(
  mode="classification",
  cost=tune(),
  degree=1
) %>%
  set_engine("kernlab")

# svm_grid<-grid_regular(cost(range=c(50,550)),
#                        levels=10,
#                        original = F)

svm_grid = tibble(cost = seq(0,70,by=2))

wkfl_svm = wkfl %>%
  add_model(genome_svm)

doParallel::registerDoParallel()
svm_start = Sys.time()
svm_res = tune_grid(wkfl_svm,
                    resamples = cv_splits,
                    grid = svm_grid,
                    metrics = metric_set(accuracy),
                    control = control_grid(save_pred = T))

svm_end = Sys.time()

svm_metrics = collect_metrics(svm_res)

ggplot(data = svm_metrics,
       aes( x = cost, y = mean)) +
       geom_line()

#get best svm
best_tuning = svm_res %>% 
  select_best(metric = "accuracy")

wkfl_best_svm = finalize_workflow( wkfl_svm , best_tuning)

final_svm = fit (wkfl_best_svm , data = genome_train)

# KNN model ----

genome_knn = nearest_neighbor(
  mode = "classification",
  neighbors= tune(),
  weight_func = "epanechnikov", 
  dist_power = 1
) %>%
  set_engine("kknn")

knn_grid = tibble(neighbors= seq(1, 151 , by = 5))

wkfl_knn = wkfl %>%
  add_model(genome_knn)

doParallel::registerDoParallel()
knn_start=Sys.time()

knn_res = tune_grid(
  wkfl_knn,
  resamples = cv_splits, 
  grid = knn_grid, 
  metrics = metric_set(accuracy),
  control = control_grid(save_pred = TRUE)
)
knn_fin = Sys.time()

knn_metrics = collect_metrics(knn_res)

ggplot( data = knn_metrics , 
        aes( x = neighbors, y = mean)) +
  geom_point()


# Model Assessment

#obtain soft classification on testing set
rf_pred<-predict(final_rf,genome_test,type="prob")
lr_pred<-predict(final_lr,genome_test,type="prob")
svm_pred<-predict(final_svm,genome_test,type="prob")
model_preds = list(lr_pred, rf_pred, svm_pred )

model_names<-c("logistic regression","random forest","support vector machine")

pacman::p_load(cowplot)
#grab true responses in test set
truth<-genome_test$sweep

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
