#script for making parallel coordinates plots for each summary stat

#load in data
genomes = read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/split_snp/0.25ds_snp_set.csv")

#downsample data to make plots clearer
set.seed(1707)
genomes = dplyr::sample_n(genomes, nrow(genomes)*0.1)

genome_SS <- genomes %>%
  dplyr::select(ID, sweep,H_1:h123_11,severity)

taj_D = rep(NA,11)
for(i in 1:11){
  taj_D[i]= paste("TajD_",i,sep="")
}

colnames(genome_SS) [which(colnames(genome_SS)=="D_1"):which(colnames(genome_SS)=="D_11")] <- taj_D 

#convert to long form
genome_SS_long = 
  genome_SS %>%
  pivot_longer(H_1:h123_11)

labs  <- glue::glue("TajD[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "TajD")) %>% 
  mutate(name = factor(name, levels = str_c("TajD_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.005) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) +
  facet_wrap("severity")

genome_SS_long %>% 
  filter(str_detect(name, "TajD")) %>% 
  mutate(name = factor(name, levels = str_c("TajD_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) 

labs  <- glue::glue("H[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "H")) %>% 
  mutate(name = factor(name, levels = str_c("H_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h1[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h1_")) %>% 
  mutate(name = factor(name, levels = str_c("h1_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h2[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h2")) %>% 
  mutate(name = factor(name, levels = str_c("h2_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h12[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h12_")) %>% 
  mutate(name = factor(name, levels = str_c("h12_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("h123[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "h123")) %>% 
  mutate(name = factor(name, levels = str_c("h123_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))

labs  <- glue::glue("w_max[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "w_max")) %>% 
  mutate(name = factor(name, levels = str_c("w_max_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.015) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) +
  scale_y_continuous(trans='log10')

labs  <- glue::glue("Zns[{1:11}]")
genome_SS_long %>% 
  filter(str_detect(name, "Zns")) %>% 
  mutate(name = factor(name, levels = str_c("Zns_", 1:11))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  geom_smooth(aes(group = sweep), se = FALSE) + 
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs))
