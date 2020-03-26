pacman::p_load(tidyverse)
library(tidymodels)

#exploratory data analysis

df<-read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/bt_cpop.csv")
df$sweep = as.factor(df$sweep)
df$s_coef=as.factor(df$s_coef)

#Data cleaning. Take out bottleneck info for now. 
genome_SS  <- df %>% 
  select(sweep, H_1:h123_11)
genome_SS

#PCA

#Using all predictors for PCs

genome_PCA_all  <- 
  df %>% 
  select(sweep, H_1:h123_11) %>% 
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

gpca <- prcomp(genome_SS[,-1], scale=T)
factoextra::fviz_eig(gpca)

# Using just the central 6 windows

genome_PCA_6  <- 
  df %>% 
  select(sweep, H_6, D_6, h1_6) %>% 
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

#Parallel Coordinates Plot

genome_SS <- df %>%
  select(ID, sweep,H_1:h123_11)

#convert to long form
genome_SS_long = 
  genome_SS %>%
  pivot_longer(H_1:h123_11)

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

labs  <- glue::glue("H[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "H")) %>% 
  mutate(name = factor(name, levels = str_c("H_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h1[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h1_")) %>% 
  mutate(name = factor(name, levels = str_c("h1_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h2[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h2")) %>% 
  mutate(name = factor(name, levels = str_c("h2_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))



