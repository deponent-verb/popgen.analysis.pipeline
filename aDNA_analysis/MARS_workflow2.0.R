pacman::p_load(tidyverse,vip,drake,earth)
library(tidymodels)

#read in data
#ancient_genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
ancient_genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
# ancient_genomes$sweep <- ifelse(ancient_genomes$sweep=="hard",1,0)
ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)


#make sure hard is the success level, neutral is reference level (R takes first level as ref.)
# ancient_genomes$sweep <- as.factor(ancient_genomes$sweep)
# ancient_genomes$sweep = relevel(ancient_genomes$sweep, "neutral")
# contrasts(ancient_genomes$sweep)

#read in functions. Must not be done outside because drake calls functions from
#external environment. 
source("~/Documents/GitHub/popgen.analysis.pipeline/aDNA_analysis/auc_aDNA.R")
source("~/Documents/GitHub/popgen.analysis.pipeline/aDNA_analysis/MARS_fit.R")

genomes_group <- function(){
  #go make a new factor of combined factor variables
  ancient_genomes %>%
    mutate(tech = interaction(impute_method,denoise_method)) %>%
    split(f = .$tech)
}

fit_model <- function(genomes_group){
  genomes = genomes_group[[1]]
  
  set.seed(12)
  #split data into training and test set
  genome_split <- initial_split(genomes,prop=0.8,strata = missing_rate)
  genome_train = training(genome_split)
  genome_test = testing(genome_split)
  
  #fit MARS model using specified ml workflow
  final_workflow = MARS_fit(genomes)
  
  #compute AUC for each set of missing rates and s_coef
  tech = genomes$tech %>% unique()
  auc_res = auc_aDNA(fitted_model = final_workflow, 
                     test_data = genome_test,
                     recipe = std_recipe) %>%
    mutate(tech = tech)
}

plan <- drake_plan(
  genomes = genomes_group(),
  model = target(fit_model(genomes), dynamic = map(genomes))
)

make(plan)
readd(model)

results = readd(model)

#generate AUC tables

table = results %>%
  split(f =.$tech)

sapply(table, function(x){mean(x$.estimate)})

#generate AUC plots

ggplot(data = results, aes(x = missing_rate, y = .estimate, color = factor(s_coef))) +
  geom_point() +
  geom_line(aes(group = factor(s_coef))) +
  ylab("AUC") +
  facet_wrap("tech", scales = "free_y") +
  scale_color_discrete(name = "technique")

results_zero = results %>% 
  filter( (tech == "zero.cluster") | (tech == "zero.fixed_cluster") | (tech =="zero.none"))

ggplot(data = results_zero, aes(x = missing_rate, y = .estimate, color = tech)) +
  geom_point() +
  geom_line(aes(group = tech)) +
  ylab("AUC") +
  facet_wrap(~s_coef)+
  scale_color_discrete(name = "technique")

results_random = results %>% 
  filter( (tech == "random.cluster") | (tech == "random.fixed_cluster") | (tech =="random.none"))

ggplot(data = results_random, aes(x = missing_rate, y = .estimate, color = tech)) +
  geom_point() +
  geom_line(aes(group = tech)) +
  ylab("AUC") +
  facet_wrap(~s_coef)+
  scale_color_discrete(name = "technique")

