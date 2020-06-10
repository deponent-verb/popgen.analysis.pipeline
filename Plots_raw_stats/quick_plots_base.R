library(ggplot2)
library(data.table)
library(cowplot)
library(tidymodels)
library(tidyverse)
library(skimr)
library(recipes)
library(forcats)

setwd("~/repos/popgen.analysis.pipeline/Plots_raw_stats/")


# Looking at raw data

raw = fread(input = "../data/notrim_base.csv")

raw[, t1 := bottle_time1]
raw[, s := bottle_size1]
raw[, t2 := bottle_time2]
raw[, Duration := t1-t2]
raw[, Onset := t2]
raw[, Strength := s]

# Melting the data

raw_m = melt(raw, measure.vars = 8:129, variable.name = "Statistic_Window", value.name = "Value")
raw_m[, Window := as.numeric(gsub(Statistic_Window, pattern = ".+_(\\d+)$", replacement = "\\1"))]
raw_m[, Statistic := gsub(Statistic_Window, pattern = "^(.+)_\\d+$", replacement = "\\1")]



# Points and smooth
# Const size
p1 = ggplot(raw_m[severity == 0], aes(Window, Value, col=factor(s_coef))) + geom_point(alpha=0.1) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

# Difficult bottlenecks:
p2 = ggplot(raw_m[demography == "t1:1680_s:0.05_t2:1600"], aes(Window, Value, col=factor(s_coef))) + geom_point(alpha=0.1) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"))


# The same plots without the points
# Const size
p1 = ggplot(raw_m[severity == 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

# Difficult bottleneck:
p2 = ggplot(raw_m[severity == 80], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_wrap("Statistic", scales = "free_y", ncol=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"))


# Plot with violins:
# Const size
p1 = ggplot(raw_m[severity == 0], aes(factor(Window), Value, fill=factor(s_coef))) + geom_violin(scale = "width") + facet_wrap("Statistic", scales = "free_y", nrow=2)

# Difficult bottleneck:
p2 = ggplot(raw_m[severity == 80], aes(factor(Window), Value, fill=factor(s_coef))) + geom_violin(scale = "width") + facet_wrap("Statistic", scales = "free_y", nrow=2)

plot_grid(p1, p2, labels = c("Const", "Bottleneck"), nrow = 2)





# Look at stats for all bottlenecks, central window 6, violin plots

# Just TajD
ggplot(raw_m[Statistic == "D" & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2)

# Just Fay&Wu's H
ggplot(raw_m[Statistic == "H" & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2)

# Number SNPs
ggplot(raw_m[Statistic == "block_snp_length" & severity != 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_grid(Onset~Duration + Strength)

ggplot(raw_m[Statistic == "base_length" & severity != 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_grid(Onset~Duration + Strength) + geom_hline(yintercept = 1/11, lty=2)


# All stats

listOfStats = raw_m[,unique(Statistic)]

p_list = lapply(listOfStats, function (x) {
  ggplot(raw_m[Statistic == x & Window == 6 & severity != 0], aes(factor(s_coef), Value, fill=factor(Strength))) + geom_violin(scale = "width") + facet_grid(Onset~Duration) + geom_hline(yintercept = 0, lty=2) + ggtitle(label=x)
})

plot_grid(plotlist = p_list, nrow = 4)


# Look at stats for all bottlenecks, all windows, smooth interpolation

p_list = lapply(listOfStats, function (x) {
  ggplot(raw_m[Statistic == x & severity != 0], aes(Window, Value, col=factor(s_coef))) + geom_smooth() + facet_grid(Onset~Duration + Strength) + ggtitle(label=x)
})

plot_grid(plotlist = p_list, nrow = 2)













######################
######################
######################
#### PCA

genomes = fread(input = "../data/0.25ds_set.csv")

#Data cleaning. Take out bottleneck info and selection coefficient for now. ----
genome_SS  <- genomes %>% 
  #filter (demography == 'cpop') %>%
  filter (s_coef %in% c(0, 0.005)) %>%
  dplyr::select(sweep, H_1:h123_11)
genome_SS

skimr::skim(genome_SS)

 #PCA ----

#We will make a recipe to transform our data into their PCs.

pc_all <- recipe(data = genome_SS, sweep + demography ~.) %>% #set our group to be sweep
  step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
  step_pca(all_predictors()) %>% #compute PCs
  prep() #give it the data to compute PC's.

tidy_pca <- tidy(pc_all, 2) #take the loadings on each PC. Entered 2 because our recipe had 2 steps.

tidy_pca = data.table(tidy_pca)
tidy_pca[, Window := as.numeric(gsub(terms, pattern = ".*_(\\d+)$", replacement = "\\1"))]
tidy_pca[, Stat := gsub(terms, pattern = "(.*)_\\d+$", replacement = "\\1")]

ggplot(tidy_pca[tidy_pca$component %in% c("PC1","PC2"),], aes(Window, value, col=Stat)) + geom_line() + geom_point() +facet_grid(component ~ Stat) + geom_hline(yintercept = 0, lty=2)




#We can now transform our data into PCs, using our recipe. 

PCA_all_demos = list()
PCA_loadings_all_demos = list()
for (demo in unique(genomes$demography)) {
  genome_SS  <- genomes %>% 
    filter (demography == demo) %>%
    #filter (s_coef %in% c(0, 0.005)) %>%
    dplyr::select(demography, sweep, H_1:h123_11)
  pc_all <- recipe(data = genome_SS, sweep + demography ~.) %>% #set our group to be sweep
    step_normalize(all_predictors()) %>%  #standardise all predictors to have mean 0, var 1.
    step_pca(all_predictors()) %>% #compute PCs
    prep() #give it the data to compute PC's.
  
  tidy_pca <- data.table(tidy(pc_all, 2)) #take the loadings on each PC. Entered 2 because our recipe had 2 steps.
  tidy_pca[, Window := as.numeric(gsub(terms, pattern = ".*_(\\d+)$", replacement = "\\1"))]
  tidy_pca[, Stat := gsub(terms, pattern = "(.*)_\\d+$", replacement = "\\1")]
  
  if (tidy_pca[Stat == "Zns" & component == "PC1", value][1] < 0) tidy_pca[, value := -value]   # Switch sign of value of loadings when Zns is negative
  
  PCA_temp = bake(pc_all, genome_SS)
  
  if (demo == "cpop") {
    PCA_temp$t1 = 0
    PCA_temp$t2 = 0
    PCA_temp$s = 1
    PCA_temp$duration = 0
  } else {
    PCA_temp$t1 = as.numeric(gsub(PCA_temp$demography, pattern = "t1:(\\d+)_.*", replacement = "\\1"))
    PCA_temp$t2 = as.numeric(gsub(PCA_temp$demography, pattern = ".+t2:(\\d+)$", replacement = "\\1"))
    PCA_temp$s = as.numeric(gsub(PCA_temp$demography, pattern = ".*s:(\\d\\.\\d+)_.*", replacement = "\\1"))
    PCA_temp$duration = PCA_temp$t1 - PCA_temp$t2
  }
  tidy_pca[, c("demography", "t1", "t2", "s", "duration") := .(unique(PCA_temp$demography), unique(PCA_temp$t1), unique(PCA_temp$t2), unique(PCA_temp$s), unique(PCA_temp$duration))]
  PCA_all_demos[[demo]] = PCA_temp
  PCA_loadings_all_demos[[demo]] = tidy_pca
}


genome_PCA_all = rbindlist(PCA_all_demos)
genome_PCA_loadings_all = rbindlist(PCA_loadings_all_demos)

#We can now plot our data with the first 2 PC's. Group by sweep.

genome_PCA_all %>% 
  filter (demography != "cpop") %>% 
  #specify PC's as axis. Group by sweep.
  ggplot(aes(x = PC1, y = PC2, col = sweep)) + 
  geom_point(alpha = 0.01) + 
  geom_density_2d() + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("PCA plot using all predictors") +
  facet_grid(t2 ~ duration + s) +
  theme(plot.title = element_text(hjust = 0.5)) 

# PC1 loadings
genome_PCA_loadings_all %>% 
  filter (duration != 0 & component %in% c("PC1")) %>% 
  #specify PC's as axis. Group by sweep.
  ggplot(aes(x = Window, y = value, col = Stat)) + 
  geom_line() + 
  geom_point() +
  #scale_color_brewer(palette = "Set2") +
  ggtitle("PCA plot using all predictors") +
  facet_grid(t2 + duration + s ~ Stat) +
  geom_hline(yintercept = 0, lty=2) +
  theme_gray()

# PC2 loadings
genome_PCA_loadings_all %>% 
  filter (duration != 0 & component %in% c("PC2")) %>% 
  #specify PC's as axis. Group by sweep.
  ggplot(aes(x = Window, y = value, col = Stat)) + 
  geom_line() + 
  geom_point() +
  #scale_color_brewer(palette = "Set2") +
  ggtitle("PCA plot using all predictors") +
  facet_grid(t2 + duration + s ~ Stat) +
  geom_hline(yintercept = 0, lty=2) +
  theme_gray()

