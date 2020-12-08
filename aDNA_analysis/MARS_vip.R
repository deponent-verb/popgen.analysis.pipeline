#script to investigate VIP for MARS model

pacman::p_load(tidyverse,vip,drake,earth)
library(tidymodels)

#read in data
#ancient_genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)


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
  
  final_workflow %>%
      pull_workflow_fit() %>%
      vip(num_features = 10, method = "firm", train = genomes) +
      ggtitle(tech)
  
  # final_workflow %>%
  #   pull_workflow_fit() %>%
  #   vip(num_features = 10, method = "firm", train = genomes) %>%
  #   mutate(tech = tech)
}

plan <- drake_plan(
  genomes = genomes_group(),
  model = target(fit_model(genomes), dynamic = map(genomes))
)

make(plan)
readd(model)

#check pdp plots for preferred method

genomes = ancient_genomes %>%
  mutate(tech = interaction(impute_method,denoise_method)) %>%
  dplyr::filter(tech == "random.fixed_cluster")

final_workflow = MARS_fit(genomes)
