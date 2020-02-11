#load packages
pacman::p_load("popgen.tools","tidyverse","ggplot2","GGally","caret")

#read in data

hard<-readRDS("~/work/MPhil/data/hard.rds")
hard2<-readRDS("~/work/MPhil/data/hard2.rds")
neutral<-readRDS("~/work/MPhil/data/neutral.rds")
neutral2<-readRDS("~/work/MPhil/data/neutral2.rds")
data<-c(hard,hard2,neutral,neutral2)
saveRDS(data,file = "~/work/MPhil/data/toy_data.rds")
# data<-readRDS("~/work/MPhil/data/toy_data.rds")

#snp distribution----

#check SNP distribution
snp_dist<-snp_count(data)

#check snp distribution using boxplots
ggplot(snp_dist,aes(sweep_type,SNP))+geom_boxplot()

ggplot(data=snp_dist, aes(x=SNP, color=sweep_type))+ geom_density()

#bimodal distribution due to two different selection coefficients. 

temp<-snp_dist %>% filter(sweep_type=="hard") %>% filter (s==0.01) %>% select(SNP)
low_mean<-mean(temp$SNP) %>% round()
low_std<-sd(temp$SNP) %>% round()
snp_cutoff<-low_mean-2*low_std #2600

#generate the dataframe ----

#data<-c(hard,soft)
#saveRDS(data,file="~/work/MPhil/data/toy_data.rds")

data<-readRDS("~/work/MPhil/data/toy_data.rds")
df<-generate_df(sim_list = data,win_split = 11,snp=snp_cutoff,form="wide")

#check if there's any NAs. That would make me sad. 
apply(df, 2, function(x) any(is.na(x)))

readr::write_csv(df,path="./data/toy_df.csv")

## Read in dataframe with raw data

df<-read_csv("./data/toy_df.csv")

df<-as_tibble(df)
df$sweep<-df$sweep %>% as.factor()
str(df)

#drop ID column
df<-subset(df,select=-ID)

#there are no NAs!
#df<-df %>% drop_na()

# exploratory data analysis----

#parallel coords plot

#create df to store means of each variable for both classes
mean_values<- df %>% group_by(sweep) %>% 
  summarise_all(mean) %>%
  pivot_longer(-sweep,names_to = "variable", values_to = "value")

#check mean computations
df %>% dplyr::filter(sweep=="hard") %>% summarise_all(mean)
df %>% dplyr::filter(sweep=="neutral") %>% summarise_all(mean)

#grabbing the right values for each variable (not complete)

p<-ggparcoord(data=df,columns = (3):(13),groupColumn="sweep",scale="globalminmax",alphaLines = 0.1) +
  geom_point(data=mean_values[2:12,],aes(x=variable, y=value, color=sweep),size=3,inherit.aes = F) +
  geom_point(data=mean_values[69:79,],aes(x=variable, y=value, color=sweep),size=3,inherit.aes = F)

p2<-ggparcoord(data=df,columns = (14):(24),groupColumn="sweep",scale="globalminmax",alphaLines = 0.1) +
  geom_point(data=mean_values[13:23,],aes(x=variable, y=value, color=sweep),size=3,inherit.aes = F) +
  geom_point(data=mean_values[80:90,],aes(x=variable, y=value, color=sweep),size=3,inherit.aes = F)

p3<-ggparcoord(data=df,columns = (26):(36),groupColumn=2,scale="globalminmax")
p4<-ggparcoord(data=df,columns = (37):(47),groupColumn=2,scale="globalminmax")
p5<-ggparcoord(data=df,columns = (48):(58),groupColumn=2,scale="globalminmax")
p6<-ggparcoord(data=df,columns = (59):(69),groupColumn=2,scale="globalminmax")

#remove ID column
df<-select(df,-c(ID))

#PCA 
library("ggfortify")
df_nores<-select(df,-c(sweep,s_coef))
autoplot(prcomp(df_nores,scale=T))

autoplot(prcomp(df_nores,scale=T),data=df,colour='sweep')


#check for near 0 variance predictors
pacman::p_load(caret)

nzv<-nearZeroVar(df,saveMetrics = TRUE)
#there are no problematic variables, not counting ID. We want F on both zeroVar and nzv. 

#any problematic variables can be removed as such
#filter_df<-df[,-nzv]

#convert data into a design matrix
df_m<-dummyVars(sweep~.,data=df)
head(predict(df_m, newdata = df))


#identify correlated predictors
df_pred<-select(df,-c(sweep,s_coef))
des_cor<-cor(df_pred)
high_corr<-findCorrelation(des_cor,cutoff = 0.75)
filtered_pred<-df_pred[,-high_corr]

#check that we removed predictors with corr >0.75
des_cor2<-cor(filtered_pred)
summary(des_cor2[upper.tri(des_cor2)])

#check for linear combinations
comboInfo<-findLinearCombos(filtered_pred)
filtered_pred<-filtered_pred[,-comboInfo$remove]

## scaling ----

#substract mean and divide by std for each column. Note that we probably want to do this over a row/window. 
preProcValues<-preProcess(filtered_pred,method=c("center","scale"))

#make new transformed set of training and test data
dataTransformed<-predict(preProcValues,filtered_pred)
final_tdata<-cbind(df$sweep,dataTransformed)
colnames(final_tdata)[1]<-"sweep"

#trainTransformed<-predict(preProcValues,train_data)
#testTransformed<-predict(preProcValues,test_data)

#partition dataset

set.seed(1688)
train_size=900
sample_size=nrow(final_tdata)
inTrain<-sample(sample_size,train_size)

train_data<-final_tdata[inTrain,]
test_data<-final_tdata[-inTrain,]

# model fitting ----

#do 10 fold CV, 3 times for each model
train.control<-trainControl(method="repeatedcv", number=10, repeats=3,classProbs = TRUE)

#logistical regression, boosted
grid<-expand.grid(C=seq(5,20,by=5))
lrfit<-train(sweep~., data=train_data, method="LogitBoost",
             nIter=grid,metric="Accuracy",trControl=train.control )

lr.preds<-predict(lrfit, test_data) 
confusionMatrix(lr.preds, test_data$sweep)

#xgboost

tune.grid <- expand.grid(eta = c(0.05, 0.075, 0.1), nrounds = c(50, 75, 100),
                         max_depth = 6,
                         min_child_weight = c(2.0, 2.25, 2.5), colsample_bytree = c(0.3, 0.4, 0.5), gamma = 0,
                         subsample = 1)

xg.fit <- train(sweep ~ .,
                data = train_data,
                method = "xgbTree", tuneGrid = tune.grid, trControl = train.control)


xg.preds<-predict(xg.fit,test_data) 
confusionMatrix(xg.preds, test_data$sweep)

#SVM (Linear)

grid<-expand.grid(C=seq(0.5,1.5,by=0.2))

svm.fit<-train(sweep~., data=train_data,method="svmLinear",
               tuneGrid=grid,metric="Accuracy",trControl=train.control)

svm.preds<-predict(svm.fit,test_data,type="prob")
confusionMatrix(svm.preds,test_data$sweep)

#SVM (poly)

grid<-expand.grid(degree=seq(1,4,by=1),scale=seq(0.1,1,by=0.2),C=seq(0.5,1.5,by=0.2))

psvm.fit<-train(sweep~., data=train_data,method="svmPoly",
               tuneGrid=grid,metric="Accuracy",trControl=train.control)

psvm.preds<-predict(svm.fit,test_data)
confusionMatrix(psvm.preds,test_data$sweep)

#Random Forest

grid<-expand.grid(mtry=seq(15,35,by=5))

rf.fit<-train(sweep~.,data=train_data,method="parRF",tuneGrid=grid,metric="Accuracy",trControl=train.control)
rf.preds<-predict(rf.fit,test_data,type = "prob")
confusionMatrix(rf.preds,test_data$sweep)

#Model Comparison
pacman::p_load(yardstick)

truth<-test_data$sweep
rf_roc<-cbind(truth,rf.preds)
roc_curve(rf_roc,truth=truth,estimate=hard)

roc_curve(rf_roc,truth=truth,estimate=hard) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity)) + 
  geom_path() +
  geom_abline(lty=3) +
  coord_equal() +
  theme_bw()

svm_roc<-cbind(truth,svm.preds)

roc_curve(svm_roc,truth=truth,estimate=hard) %>% 
  ggplot(aes(x = 1 - specificity, y = sensitivity)) + 
  geom_path() +
  geom_abline(lty=3) +
  coord_equal() +
  theme_bw() 

rf_out<-roc_curve(rf_roc,truth=truth,estimate=hard) %>% 
  mutate(model="rf")

svm_out<-roc_curve(svm_roc,truth=truth,estimate=hard) %>%
  mutate(model="svm")

models_out<-bind_rows(rf_out,svm_out)

models_out %>% 
  #get individual roc curves for each model
  group_by(model) %>% 
  ggplot(
    aes(x=1-specificity,y=sensitivity,color=model)
    ) +
  geom_line(size=1.1) +
  geom_abline(slope=1, intercept = 0, size=0.2) +
  #fix aspect ratio to 1:1 for visualization
  coord_fixed() +
  theme_cowplot()

#AUC
  

#stop here
#yardstick for ROC curve,

#comparison, 

#compute time
#ROC
#explanatory
#confusion matrix (standard cutoff)
#variable of importance
#keras DL

################################################

#sanity checking of summary statistics (raw data) ----

#index of starting column with summary stats
sum_start=4

#number of splits in a window
wins=11

ggparcoord(data=df, columns=15:25,groupColumn = 1)

#unscaled parcoords
for(i in 1:5){
  p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
  print(p)
}

i=1
 p<-ggparcoord(data=df%>%filter(sweep=="neutral"),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df%>%filter(sweep=="soft"),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df%>%filter(sweep=="hard")%>%filter(s_coef==0.01),columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)
 print(p)
 
 temp<-df%>%filter(sweep=="hard") %>% head(50)
 
 p<-ggparcoord(data=temp[,15:25],columns = 1:11)
 plot(p)
 
 temp<-df %>% filter(sweep=="hard") 
 
 ggplot(data=temp[,15:25],aes())+geom_line()

#scaled parcoords
for(i in 1:5){
  p<-ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1,scale="center")
  print(p)
}

#check boxplots

ggplot(data=df,aes(x =sweep,y=D_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h1_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h2_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h12_6))+geom_boxplot()+scale_y_log10()
ggplot(data=df,aes(x =sweep,y=h123_6))+geom_boxplot()+scale_y_log10()



#### outdated code
i<-rep(0:5)
ggparcoord(data=df,columns = (sum_start+i*wins):(sum_start+(i+1)*wins-1),groupColumn=1)

df %>% pivot_longer(H1:h123_10) %>% 
  mutate(batch = str_remove_all(name, '\\d')) %>% 
  ggplot(aes(name, value)) + geom_line(aes(group = ID)) + 
  facet_wrap(~batch, scales = "free")

ggplot(data=df, aes(x=dist,y=D)) + geom_line(aes(color=sweep)) + scale_y_continuous(trans="log10")


DF<-readRDS("~/work/MPhil/data/df.rds")

sim<-data[[414]]
sub_win(sim$genomes,2)


#generate new IDs so that they are not informative. 

size<-nrow(df)/10
x<-rep(1:size)
ID<-sample(x,size=size)

foo<-function(y,it){
  ans<-rep(y,it)
  return(ans)
}


new_ID<-map2(ID,10,foo)
test<-unlist(new_ID)

#The interpretibility of the model is helpful for understanding the underlying biology that's giving you that answer.
#Help work out what signatures are important for finding sweeps in genome. 

