#PCA for aDNA data
library(tidyverse)
library(tidymodels)
library(drake)

ancient_genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
ancient_genomes$denoise_method <- ifelse(ancient_genomes$denoise_method == "cluster", "silhouette_cluster", ancient_genomes$denoise_method)

#drake pipeline for dynamic branching
genomes_group <- function(){
  #go make a new factor of combined factor variables
  ancient_genomes %>%
    mutate(tech = interaction(impute_method,denoise_method)) %>%
    split(f = .$tech)
}

tech_pca <- function(genomes_group){
  genomes = genomes_group[[1]]
  
  genome_SS  <- genomes %>% 
    dplyr::select(s_coef, D_1:h123_5)
  
  genome_SS$s_coef <- as.factor(genome_SS$s_coef)
  
  pc_all <- recipe(data = genome_SS, s_coef ~.) %>% #set our group to be sweep
    step_naomit(all_predictors()) %>%
    step_nzv(all_predictors()) %>%
    step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
    step_pca(all_predictors()) %>% #compute PCs
    prep()
  
  tech = genomes$tech %>% unique()
  genome_PCA_all = bake(pc_all, genome_SS) %>%
    mutate(tech = tech)
}

plan <- drake_plan(
  genomes = genomes_group(),
  model = target(tech_pca(genomes), dynamic = map(genomes))
)

make(plan)
readd(model)

pca_plots = readd(model)

pca_plots %>% 
  #specify PC's as axis. Group by s_coef
  ggplot(aes(x = PC1, y = PC2, col = s_coef)) + 
  geom_point(alpha = 0.05) + 
  #geom_density_2d() + 
  scale_color_brewer(palette = "Set1") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap("tech") +
  guides(color = guide_legend(override.aes = list(alpha = 1)))


#end of dynamic branching

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

genome_SS  <- ancient_genomes %>% 
  filter (impute_method == "random") %>%
  filter (denoise_method == "none") %>%
  filter (s_coef == 0 | s_coef == 0.0125) %>%
  dplyr::select(s_coef, D_1:H_5)

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
