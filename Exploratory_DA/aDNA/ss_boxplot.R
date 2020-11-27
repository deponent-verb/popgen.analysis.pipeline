#boxplots to show how missing rate affects data separation.

pacman::p_load(tidyverse)

#genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
genomes$sweep <- as.factor(genomes$sweep)

genome_SS_long = genomes %>%
  pivot_longer(D_1:h123_5)

#SFS stats

genome_SS_long %>% 
  filter(str_detect(name, "H_")) %>%
  ggplot (aes( x = name, y = value, fill = sweep )) +
  geom_boxplot() + 
  facet_wrap("missing_rate") + 
  labs(subtitle = "title")

stats = c("D_3","H_3")
sfs_plots = list()

for(s in 1:length(stats)){
  sfs_plots[[s]] = genome_SS_long %>% 
    filter(str_detect(name, stats[s])) %>%
    ggplot (aes( x = factor(missing_rate), y = value, fill = sweep )) +
    geom_boxplot() + 
    facet_wrap("impute_method") +
    labs(x = "missing rate", title = stats[s])
}


#haplotype stats

genome_SS_long %>% 
  filter(str_detect(name, "h2_")) %>%
  ggplot (aes( x = name, y = value, fill = sweep )) +
  geom_boxplot() + 
  facet_wrap("missing_rate") + 
  labs(subtitle = "full data set")

haplo_stats = c("h1_3", "h2_3", "h12_3", "h123_3")
hplots = list()

for(h in 1:length(haplo_stats)){
  hplots[[h]] = genome_SS_long %>% 
    filter(str_detect(name, haplo_stats[h])) %>%
    mutate(tech = interaction(impute_method,denoise_method)) %>%
    ggplot (aes( x = factor(missing_rate), y = value, fill = sweep )) +
    geom_boxplot() + 
    facet_wrap("tech") +
    labs(subtitle = haplo_stats[h], x = "missing rate")
}

genome_SS_long %>% 
  filter(str_detect(name, "h2_3")) %>%
  mutate(tech = interaction(impute_method,denoise_method)) %>%
  ggplot (aes( x = factor(missing_rate), y = value, fill = sweep )) +
  geom_boxplot() + 
  facet_wrap("tech") +
  labs(subtitle = "h2_3", x = "missing rate")


#model code
DF2 <- data.frame(
  x = c(c(A1, A2, A3), c(B1, B2, B3)),
  y = rep(c("A", "B"), each = 15),
  z = rep(rep(1:3, each=5), 2),
  stringsAsFactors = FALSE)

cols <- rainbow(3, s = 0.5)

ggplot(DF2, aes(y, x, fill=factor(z))) +
  geom_boxplot()

