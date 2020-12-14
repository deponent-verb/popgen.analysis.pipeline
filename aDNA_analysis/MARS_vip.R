#script to investigate VIP for MARS model

pacman::p_load(tidyverse,vip,drake,earth)
library(tidymodels)

#read in data
#ancient_genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
ancient_genomes$sweep <- ifelse(ancient_genomes$sweep=="hard",1,0)
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)
ancient_genomes$denoise_method <- ifelse(ancient_genomes$denoise_method == "cluster", "silhouette_cluster", ancient_genomes$denoise_method)


#read in functions. Must not be done outside because drake calls functions from
#external environment. 
source("~/Documents/GitHub/popgen.analysis.pipeline/aDNA_analysis/MARS_fit.R")

genomes_group <- function(){
  #go make a new factor of combined factor variables
  ancient_genomes %>%
    dplyr::filter(impute_method == "random") %>%
    mutate(tech = interaction(impute_method,denoise_method)) %>%
    split(f = .$tech)
}

fit_model <- function(genomes_group){
  genomes = genomes_group[[1]]
  
  #fit MARS model using specified ml workflow
  final_workflow = MARS_fit(genomes)
  
  tech = genomes$tech %>% unique()
  
  wkfl_rec = final_workflow %>%
    pull_workflow_prepped_recipe()
  
  final_workflow %>%
    pull_workflow_fit() %>%
    #need to supply transformed training data for vip
    vip(method = "firm", train = bake(wkfl_rec,genomes)) +
    ggtitle(tech)
  
  # final_workflow %>%
  #   pull_workflow_fit() %>%
  #   vip(num_features = 10, method = "firm", train = genomes) %>%
  #   mutate(tech = tech)
  
  # final_workflow %>%
  #     pull_workflow_fit() %>%
  #     vip(num_features = 10, method = "firm", train = genomes) +
  #     ggtitle(tech)
}

plan <- drake_plan(
  genomes = genomes_group(),
  model = target(fit_model(genomes), dynamic = map(genomes))
)

make(plan)
plots = readd(model)

par(mfrow=c(1,3))
plots[[1]]
plots[[2]]

grid.arrange(plots[[1]],plots[[2]],plots[[3]], ncol = 3)
