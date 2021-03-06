#parallel coordinates plots for aDNA data
pacman::p_load(tidyverse)

#genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
genomes = read_csv("~/Documents/GitHub/popgen.analysis.pipeline/data/cleaned_aDNA.csv")
genomes$sweep <- as.factor(genomes$sweep)

genome_SS_long = genomes %>%
  filter(impute_method == "random") %>%
  #filter(dmg_rate == 0.02) %>%
  filter(denoise_method == "none") %>%
  pivot_longer(D_1:h123_5)

#Plot
labs  <- glue::glue("D[{1:5}]")

genome_SS_long %>% 
  filter(str_detect(name, "D_")) %>% 
  mutate(name = factor(name, levels = str_c("D_", 1:5))) %>% 
  ggplot(aes(x = name, y = value, col = sweep)) + 
  geom_line(aes(group = ID), alpha = 0.01) + 
  #geom_smooth(aes(group = sweep), se = FALSE) + 
  stat_summary(fun=mean, aes(col=sweep, group = sweep), geom="line", size = 1) +
  #stat_summary(fun=mean, colour="red", geom="line", size = 3) # draw a mean line in the data
  scale_color_brewer(palette = "Set1") + 
  labs(x = "Summary statistic", y = "Observed value", col = "Type of sweep") + 
  scale_x_discrete(labels = parse(text = labs)) +
  facet_wrap("missing_rate")

#if you get Error in grid.Call(C_textBounds, as.graphicsAnnot(x$label), x$x, x$y,  : 
#polygon edge not found, make sure you run the labs line as well.

#rerun if you get the same error

#as missing rate increases, the haplotype statistics become less informative
