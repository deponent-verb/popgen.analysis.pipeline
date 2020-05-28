pacman::p_load(tidyverse)
library(tidymodels)

# load in cleaned data from data cleaning script
genomes = read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/0.25ds_set.csv")

#Data cleaning. Take out bottleneck info and selection coefficient for now. ----
genome_SS  <- genomes %>% 
  filter (demography == 'cpop') %>%
  dplyr::select(sweep, H_1:h123_11)
genome_SS

skimr::skim(genome_SS)

#PCA ----

#We will make a recipe to transform our data into their PCs.

pc_all <- recipe(data = genome_SS, sweep ~.) %>% #set our group to be sweep
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep() #give it the data to compute PC's.

tidy_pca <- tidy(pc_all, 2) #take the loadings on each PC. Entered 2 because our recipe had 2 steps.

#plot the PC loadings for the top 4 PC's. 
tidy_pca %>%
  filter(component %in% c("PC1","PC2")) %>% #plot first 4 PC's
  mutate(component = fct_inorder(component)) %>% #make sure PC's are in order
  ggplot (aes(value, terms, fill=terms)) + 
  geom_col(show.legend = F) + 
  facet_wrap(~component) + 
  labs(y=NULL)

#Plot the top 20 loadings by absolute value for each PC 
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

#PC loadings in depth ----
#Group the summary statistics and look at their respective loadings on the top PC's.

#Group statistics on via SFS, LD, haplotype

pca_load <- tidy_pca %>%
  #grab the first 2/3 letters of the terms to classify them 
  mutate(type = case_when(stringr::str_extract(terms, "^.{2}") %in% c("H_", "D_") ~ 'SFS',
                          stringr::str_extract(terms, "^.{2}") %in% c("h1", "h2") ~ 'Haplotype',
                          stringr::str_extract(terms, "^.{3}") %in% c("Zns", "LD_","w_m") ~ 'LD',
                          TRUE ~ 'Other'
  )) %>%
  #group via window
  mutate(window = stringr::str_sub(tidy_pca$terms, -2, -1)) %>%
  #group via statistic
  mutate(stat = case_when(stringr::str_extract(terms, "^.{2}") == "H_" ~ 'H',
                          stringr::str_extract(terms, "^.{2}") == "D_" ~ 'D',
                          stringr::str_extract(terms, "^.{4}") == "LD_a" ~ 'LD_avg',
                          stringr::str_extract(terms, "^.{4}") == "LD_m" ~ 'LD_max',
                          stringr::str_extract(terms, "^.{3}") == "Zns" ~ 'Zns',
                          stringr::str_extract(terms, "^.{2}") == "w_" ~ 'w_max',
                          stringr::str_extract(terms, "^.{3}") == "h1_" ~ 'h1',
                          stringr::str_extract(terms, "^.{3}") == "h2_" ~ 'h2',
                          stringr::str_extract(terms, "^.{4}") == "h12_" ~ 'h12',
                          stringr::str_extract(terms, "^.{4}") == "h123" ~ 'h123'
                          ))

#boxplot of loadings grouped by type
pca_load %>%
  filter(component %in% c("PC1","PC2","PC3","PC4")) %>%
  ggplot(aes(y = value, color = type)) + 
  geom_boxplot() + 
  facet_wrap(~component) + 
  ggtitle("PCA on cpop data")

#boxplot of loadings grouped by statistic
pca_load %>%
  filter(component %in% c("PC1","PC2","PC3","PC4")) %>%
  ggplot(aes(y = value, color = stat, x = stat)) + 
  geom_boxplot() + 
  facet_wrap(~component) + 
  ggtitle("PCA on cpop data")

#boxplot of loadings grouped by window
pca_load %>%
  filter(component %in% c("PC1","PC2","PC3","PC4")) %>%
  ggplot(aes(y = value, color = window)) + 
  geom_boxplot() + 
  facet_wrap(~component) + 
  ggtitle("PCA on cpop data")
  
#Ask Jono how to interpret this plot. Does a high magnitude value mean that the predictor does a good job separating the data?
#Need Jono's expert coding skills to tidy this up. (e.g. fix the horizontal axis, add color code)