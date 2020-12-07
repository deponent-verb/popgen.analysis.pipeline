#quick LDA plot
library(MASS)

genome_SS  <- ancient_genomes %>% 
  filter (impute_method == "random") %>%
  filter (denoise_method == "none") %>%
  filter (s_coef == 0 | s_coef == 0.1) %>%
  dplyr::select(sweep, D_1:H_5)

#genome_SS$s_coef <- as.factor(genome_SS$s_coef)
genome_SS$sweep <- as.factor(genome_SS$sweep)

#CV=T for jack-knifed estimates
lda_fit = lda(sweep~., data = genome_SS)
plot(lda_fit,dimen = 2,nbins = 100)

#devtools::install_github("fawda123/ggord")
library(ggord)
ggord(lda_fit, genome_SS$sweep,alpha = 0.01)
