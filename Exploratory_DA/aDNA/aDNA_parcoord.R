#parallel coordinates plots for aDNA data
pacman::p_load(tidyverse)

genomes = read_csv("./data/cleaned_aDNA_nodeam.csv")
genomes$sweep <- as.factor(genomes$sweep)

genome_SS_long = 
  genomes %>%
  pivot_longer(D_1:h123_5)

#Plot
labs  <- glue::glue("h1[{1:5}]")

genome_SS_long %>% 
  filter(str_detect(name, "h1")) %>% 
  mutate(name = factor(name, levels = str_c("h1_", 1:5))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.1) + 
  #geom_smooth(aes(group = sweep), se = FALSE) + 
  #stat_summary(fun=mean, colour="red", geom="line", size = 3) # draw a mean line in the data
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) +
  facet_wrap("missing_rate")

#if you get Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
#polygon edge not found, make sure you run the labs line as well.

#rerun if you get the same error

#as missing rate increases, the haplotype statistics become less informative
