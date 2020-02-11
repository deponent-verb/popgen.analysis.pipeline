model_names<-c("lr","xgboost","svm","rf")
model_preds<-list(lr.preds,xg.preds,svm.preds,rf.preds)

model_eval<-function(model_name,model_pred){
  data<-as.data.frame(model_pred)
  ans<-data %>% mutate(model=model_name) %>% cbind(truth)
  return(ans)
}

model_list<-list()
for(i in 1:length(model_names)){
  model_list[[i]]<-model_eval(model_names[i],model_preds[i])
}

model_out<-plyr::ldply(model_list, data.frame)
#model_list<-lapply(model_names,model_preds,FUN=model_eval)


model_out %>% 
  #get individual roc curves for each model
  group_by(model) %>% 
  roc_curve(truth=truth, estimate=hard) %>%
  ggplot(
    aes(x=1-specificity,y=sensitivity,color=model)
  ) +
  geom_line(size=1.1) +
  geom_abline(slope=1, intercept = 0, size=0.2) +
  #fix aspect ratio to 1:1 for visualization
  coord_fixed() +
  theme_cowplot()


##################################

rf_out1<-rf.preds %>% mutate(model="rf") %>% cbind(truth)
svm_out1<-svm.preds %>% mutate(model="svm")%>% cbind(truth)

model_list<-list(rf_out1,svm_out1)
temp<-plyr::ldply(model_list, data.frame)

#models_out<-bind_rows(rf_out1,svm_out1)

models_out %>% 
  #get individual roc curves for each model
  group_by(model) %>% 
  roc_curve(truth=truth, estimate=hard) %>%
  ggplot(
    aes(x=1-specificity,y=sensitivity,color=model)
  ) +
  geom_line(size=1.1) +
  geom_abline(slope=1, intercept = 0, size=0.2) +
  #fix aspect ratio to 1:1 for visualization
  coord_fixed() +
  theme_cowplot()

#AUC

a<-model_out %>%
  group_by(model) %>%
  roc_auc(truth=truth,hard,options = list(smooth = TRUE))

knitr::kable(a)

rf_out1 %>% roc_auc(truth=truth,estimate="hard",options=list(smooth=TRUE))
