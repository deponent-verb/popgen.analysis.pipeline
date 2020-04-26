pacman::p_load(tidyverse)
library(tidymodels)

#load data
snp_set<-read_csv("~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")
snp_set$sweep = as.factor(snp_set$sweep)
snp_set$s_coef=as.factor(snp_set$s_coef)
levels(snp_set$sweep) = c(1,0)

#Data cleaning. Take out bottleneck info and selection coefficient for now. ----
genome_SS  <- snp_set %>% 
  select(sweep, H_1:h123_11)
genome_SS

skimr::skim(genome_SS)

#PCA ----

#We will make a recipe to transform our data into their PCs.

pc_all <- recipe(data = genome_SS, sweep ~.) %>% #set our group to be sweep
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep() #give it the data to compute PC's.

tidy_pca <- tidy(pc_all, 2) #take the loadings on each PC. Entered 2 because our recipe had 2 steps.

tidy_pca %>%
  filter(component %in% c("PC1","PC2","PC3","PC4")) %>% #plot first 4 PC's
  mutate(component = fct_inorder(component)) %>% #make sure PC's are in order
  ggplot (aes(value, terms, fill=terms)) + 
  geom_col(show.legend = F) + 
  facet_wrap(~component) + 
  labs(y=NULL)

tidy_pca %>%
  filter(component %in% c("PC1","PC2","PC3","PC4")) %>%
  group_by(component) %>%
  top_n(20, abs(value)) %>% #find the top 10 contributors for each PC
  ungroup() %>%
  mutate(terms = drlib::reorder_within(terms, abs(value), component)) %>%
  ggplot(aes(abs(value), terms, fill = value>0)) +
  geom_col() +
  facet_wrap(~component, scales = "free_y") +
  drlib::scale_y_reordered() + 
  labs(
    x = "Absolute Value of PC Loading",
    y = NULL, fill = "Positive?"
  )


#We can now transform our data into PCs, using our recipe. 

genome_PCA_all = bake(pc_all, genome_SS)

#We can now plot our data with the first 2 PC's. Group by sweep.

genome_PCA_all %>% 
  #specify PC's as axis. Group by sweep.
  ggplot(aes(x = PC1, y = PC2, col = sweep)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("PCA plot using all predictors") +
  theme(plot.title = element_text(hjust = 0.5)) 

#get the loadings on each PC
gpca <- prcomp(genome_SS[,-1], scale=T)

#scree plot. This shows the percentage of variance explained by each PC.
factoextra::fviz_eig(gpca)

#plot loadings of each predictor on PC1
pc_load = data.frame(Value = unlist(gpca$rotation[,1]), Statistic = names(gpca$rotation[,1]))
plot(Value ~ Statistic, pc_load)

ggplot(data = pc_load, aes(x = Statistic, y = Value) ) +
  geom_point()

#Remarks. 

#Ask Jono how to interpret this plot. Does a high magnitude value mean that the predictor does a good job separating the data?
#Need Jono's expert coding skills to tidy this up. (e.g. fix the horizontal axis, add color code)