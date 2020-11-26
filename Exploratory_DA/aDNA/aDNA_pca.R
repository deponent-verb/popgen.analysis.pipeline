#PCA for aDNA data
library(tidyverse)
library(tidymodels)

genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")

#Data cleaning. Take out bottleneck info and selection coefficient for now. ----
genome_SS  <- genomes %>% 
  #filter (impute_method == "zero") %>%
  filter (denoise_method != "none") %>%
  dplyr::select(sweep, H_1:h123_5)
genome_SS

skimr::skim(genome_SS)

#PCA ----

#scree plot
rm_nz = recipe(data = genome_SS, sweep ~.) %>%
  step_nzv(all_predictors()) %>%
  prep()

scree_data = bake(rm_nz, new_data = genome_SS)
#remove response variable column
gpca <- prcomp(scree_data[,-14], scale=T)
factoextra::fviz_eig(gpca)

#pca loadings

#We will make a recipe to transform our data into their PCs.

pc_all <- recipe(data = genome_SS, sweep ~.) %>% #set our group to be sweep
  step_nzv(all_predictors()) %>% #remove non-zero variance predictors
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep() #give it the data to compute PC's.

#take the loadings on each PC. Entered 3 because our recipe had 2 steps.
tidy_pca <- tidy(pc_all, 3) 

tidy_pca %>%
  filter(component %in% c("PC1","PC2")) %>% #plot first 2 PC's
  mutate(component = fct_inorder(component)) %>% #make sure PC's are in order
  ggplot (aes(value, terms, fill=terms)) + 
  geom_col(show.legend = F) + 
  facet_wrap(~component) + 
  labs(y=NULL)

#pca plot

genome_SS  <- genomes %>% 
#  filter (impute_method == "random") %>%
  filter (denoise_method != "none") %>%
  dplyr::select(s_coef, H_1:h123_5)

genome_SS$s_coef <- as.factor(genome_SS$s_coef)

pc_all <- recipe(data = genome_SS, s_coef ~.) %>% #set our group to be sweep
  step_naomit(all_predictors()) %>%
  step_nzv(all_predictors()) %>%
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep()

genome_PCA_all = bake(pc_all, genome_SS)

genome_PCA_all %>% 
  #specify PC's as axis. Group by s_coef
  ggplot(aes(x = PC1, y = PC2, col = s_coef)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1") +
  #ggtitle("PCA plot using all predictors") +
  theme(plot.title = element_text(hjust = 0.5)) 
