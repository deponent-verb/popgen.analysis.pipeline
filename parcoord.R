pacman::p_load(tidyverse,tidymodels)

data<-read_csv("./data/toy_df.csv")
data$s_coef<-as.factor(data$s_coef)

#data cleaning
genome_SS  <- 
  data %>% 
  select(ID, sweep, D_1:h1_11)
genome_SS

# PCA 

#using D and h1 only
genome_PCA_all  <- 
  data %>% 
  select(sweep, D_1:h123_11) %>% 
  recipe(sweep ~ .) %>% 
  step_pca(all_predictors(), num_comp = 2) %>% 
  prep(training = genome_SS) %>% 
  bake(new_data = genome_SS )

genome_PCA_all %>% 
  ggplot(aes(x = PC1, y = PC2, col = sweep)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("PCA plot using all predictors") +
  theme(plot.title = element_text(hjust = 0.5))

#Just central window 6

genome_PCA_6  <- 
  data %>% 
  select(sweep, H_6, D_6, h1_6) %>% 
  recipe(sweep ~ .) %>% 
  step_pca(all_predictors(), num_comp = 2) %>% 
  prep(training = data) %>% 
  bake(new_data = data )

genome_PCA_6 %>% 
  ggplot(aes(x = PC1, y = PC2, col = sweep)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1")+
  ggtitle("PCA plot using window 6") +
  theme(plot.title = element_text(hjust = 0.5))

#comparing strong and weak selection

genome_PCA_s6  <- 
  data %>% 
  select(s_coef, H_6,D_6,h1_6) %>% 
  recipe(s_coef ~ .) %>% 
  step_pca(all_predictors(), num_comp = 2) %>% 
  prep(training = data) %>% 
  bake(new_data = data )

genome_PCA_s6 %>% 
  ggplot(aes(x = PC1, y = PC2, col = s_coef)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1")

genome_PCA_sall  <- 
  data %>% 
  select(s_coef, H_1:h123_11) %>% 
  recipe(s_coef ~ .) %>% 
  step_pca(all_predictors(), num_comp = 2) %>% 
  prep(training = data) %>% 
  bake(new_data = data )

genome_PCA_sall %>% 
  ggplot(aes(x = PC1, y = PC2, col = s_coef)) + 
  geom_point(alpha = 0.3) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1")


## Parallel Coordinates Plot

#convert to long
genome_SS_long  <- 
  genome_SS %>% 
  pivot_longer(D_1:h1_11)

labs  <- glue::glue("h1[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h1")) %>% 
  mutate(name = factor(name, levels = str_c("h1_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'

labs  <- glue::glue("D[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "D")) %>% 
  mutate(name = factor(name, levels = str_c("D_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))
## `geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'