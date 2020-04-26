#In this script, we will do some exploratory data analysis on our simulations.

#We simulated 1Mb genomes and broke them into 11 blocks. Each block has ~equal number
#of SNPs. 

#Note that for the purposes of computing the LD statistics, we did a 20% downsample 
#of each block. 

#load packages
pacman::p_load(tidyverse)
library(tidymodels)

#combining the dataframes together
neutral = read_csv("~/work/MPhil/ml_review/data/snp_split/snp_neutral_btl.csv")
hard = read_csv("~/work/MPhil/ml_review/data/snp_split/snp_hard_btl.csv")
cpop=read_csv("~/work/MPhil/ml_review/data/snp_split/snp_cpop.csv")
df = rbind(cpop,hard,neutral)
write.csv(df, "~/work/MPhil/ml_review/data/snp_split/split_snp_set.csv")

#exploratory data analysis

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
  step_center(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_scale (all_predictors()) %>%
  step_pca(all_predictors()) %>% #compute PCs
  prep(training = genome_SS) #give it the data to compute PC's.

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

#Partial Least Squares Discriminant Analysis





# Using just the central 6 windows

genome_PCA_6  <- 
  snp_set %>% 
  select(sweep, H_6, D_6,LD_avg_6, LD_max_6, w_max_6, Zns_6,
         h1_6, h2_6, h12_6, h123_6) %>% 
  recipe(sweep ~ .) %>% 
  step_pca(all_predictors(), num_comp = 2) %>% 
  prep(training = genome_SS) %>% 
  bake(new_data = genome_SS )

genome_PCA_6 %>% 
  ggplot(aes(x = PC1, y = PC2, col = sweep)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1")+
  ggtitle("PCA plot using window 6") +
  theme(plot.title = element_text(hjust = 0.5))

#PLS ----

#need plsda

gen_rec = recipe(sweep ~., data = genome_SS) %>%
  step_dummy(all_outcomes()) %>%
  step_center(starts_with("H"), starts_with("D"), starts_with("L"), 
              starts_with("w"), starts_with("Z")) %>%
  step_scale(starts_with("H"), starts_with("D"), starts_with("L"), 
             starts_with("w"), starts_with("Z")) %>%
  prep()



genome_pls <- step_pls(recipe = gen_rec, 
                       all_predictors(),
                       num_comp = 5,
                       outcome = "sweep") %>%
              prep(training = genome_SS)

# genome_PCA_all  <- 
#   snp_set %>% 
#   select(sweep, H_1:h123_11) %>% 
#   recipe(sweep ~ .) %>% 
#   step_pca(all_predictors(), num_comp = 2) %>% 
#   prep(training = genome_SS) %>% 
#   bake(new_data = genome_SS )

# genomes_recipe <- recipe(sweep ~., data=genome_train) %>%
#   step_corr(all_predictors(),threshold = 0.8) %>%
#   step_center(starts_with("H"),starts_with("D")) %>%
#   step_scale(starts_with("H"),starts_with("D")) %>%
#   prep()

pls.model = 


#Parallel Coordinates Plot

genome_SS <- snp_set %>%
  select(X1, sweep,H_1:category)

genome_SS = genome_SS %>%
  rename(
    sim_id = X1
    )

taj_D = rep(NA,11)
for(i in 1:11){
  taj_D[i]= paste("TajD_",i,sep="")
}

colnames(genome_SS) [which(colnames(genome_SS)=="D_1"):which(colnames(genome_SS)=="D_11")] <- taj_D 

#convert to long form
genome_SS_long = 
  genome_SS %>%
  pivot_longer(H_1:category)

labs  <- glue::glue("TajD[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "TajD")) %>% 
  mutate(name = factor(name, levels = str_c("TajD_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) +
  facet_wrap("category")

labs  <- glue::glue("H[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "H")) %>% 
  mutate(name = factor(name, levels = str_c("H_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h1[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h1_")) %>% 
  mutate(name = factor(name, levels = str_c("h1_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h2[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h2")) %>% 
  mutate(name = factor(name, levels = str_c("h2_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("LD_avg[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "LD_avg")) %>% 
  mutate(name = factor(name, levels = str_c("LD_avg_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("w_max[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "w_max")) %>% 
  mutate(name = factor(name, levels = str_c("w_max_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("Zns[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "Zns")) %>% 
  mutate(name = factor(name, levels = str_c("Zns_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = sim_id), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))



