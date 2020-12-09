pacman::p_load(pdp,earth,caret,tidyverse,vip)
library(tidymodels)

#pdp for our best model

#need workaround because pdp does not work with parsnip
ancient_genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
ancient_genomes$sweep <- ifelse(ancient_genomes$sweep=="hard",1,0)
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)

genomes = ancient_genomes %>%
  mutate(tech = interaction(impute_method,denoise_method)) %>%
  dplyr::filter(tech == "random.fixed_cluster")

#transform data ourselves for pdp

#Since the haplotype statistics are proportions, we don't want to normalise them.
hap_cols <- colnames(genomes)[which(colnames(genomes)=="h1_1"):which(colnames(genomes)=="h123_5")]
aDNA_dmg_cols <- c(colnames(genomes)[which(colnames(genomes)=="missing_rate"):which(colnames(genomes)=="denoise_method")],"tech")

std_recipe <- recipe(sweep ~., data = genomes) %>% #set sweep as response variable. everything else is a predictor.
  update_role(s_coef, new_role = 'demography') %>% #remove s_coef as predictor
  update_role( all_of(aDNA_dmg_cols), new_role = 'damage') %>% 
  add_role(all_of(hap_cols), new_role = 'haplotype') %>%
  step_corr(all_predictors(),threshold = 0.8) %>% #remove all highly correlated predictors
  step_normalize(all_predictors(), -has_role("haplotype")) %>% #normalize all predictors, except haplotype stats 
  prep()

baked_data = bake(std_recipe, new_data = genomes)
trans_data = baked_data[,which(colnames(baked_data)=="D_1"):which(colnames(baked_data)=="sweep")]
trans_data$tech <- NULL
trans_data$sweep <- ifelse(trans_data$sweep=="hard",1,0)


source("~/Documents/GitHub/popgen.analysis.pipeline/aDNA_analysis/MARS_fit.R")
final_workflow = MARS_fit(genomes)

final_workflow %>%
  #pull parsnip model
  pull_workflow_fit() %>%
  #pull out MARS model since pdp does not have support for parsnip
  .$fit %>%
  pdp::partial(train = trans_data, pred.var = "H_3", 
               plot = TRUE, type = "classification")

#ice plot (slow atm)
final_workflow %>%
  #pull parsnip model
  pull_workflow_fit() %>%
  #pull out MARS model since pdp does not have support for parsnip
  .$fit %>%
  pdp::partial(train = trans_data, pred.var = "H_3", 
               plot = TRUE, type = "classification", ice = TRUE)

#1:35

final_workflow %>%
  #pull parsnip model
  pull_workflow_fit() %>%
  vip()
