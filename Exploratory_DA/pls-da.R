pacman::p_load(mixOmics, tidyverse)

#load data
snp_set<-read_csv("~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")
snp_set$sweep = as.factor(snp_set$sweep)
snp_set$s_coef=as.factor(snp_set$s_coef)


#Data cleaning. Take out bottleneck info and selection coefficient for now. ----
genome_SS  <- snp_set %>% 
  dplyr::select(sweep, H_1:h123_11)
genome_SS


# To do PLS, we must separate our data into predictors and response variables
X <- genome_SS %>%
  dplyr::select(-sweep)

Y <- genome_SS$sweep

res = splsda(X,Y , scale=T)

#Plot data on the top two components
plotIndiv(res, ind.names = F , legend = T,
          ellipse = T,
          title = "PLS-plot", 
          X.label = "Comp-1", Y.label = "Comp-2")

#Display predictors that contribute >0.7 to the definition of each component.
#Further away they are from center, the higher their loading. 
plotVar(res, cutoff=0, overlap=F, cex = 2)

#plot the pls component loadings
pls_load = selectVar(res, comp =1)$value 
pls_load = setNames(cbind(rownames(pls_load), pls_load, row.names=NULL),
         c("Stat","Value"))
ggplot(data = pls_load, aes(x = Stat, y = Value)) +
  geom_point()

plotLoadings(res, comp = 1 )

#plotVar(res, cutoff=0.7, font = 0.001, overlap=F, plot = F)


#plot on classifier background
background = background.predict(res, comp.predicted = 2 , dist = "max.dist")
plotIndiv(res, comp=1:2, group = genome_SS$sweep,
          ind.names = F, title = "Max Distance",
          legend = T, background = background)

#ROC curve of pls-da classifier
auroc(res)

