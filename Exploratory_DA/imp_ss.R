#In this script, we will investigate important predictors for the LR model. 

#First we check the correlations between predictors. Certain predictors may
#be serving as proxies for others. 

pacman::p_load(corrplot,ggplot2)
genomes = read_csv("~/work/MPhil/ml_review/data/hubs_data/dataframes/split_snp/0.25ds_snp_set.csv")

#extract predictors
genomes_SS = genomes %>% 
  dplyr::select(H_1:h123_11)

corrplot(cor(genomes_SS), method="pie") #too big to see. Ask Jono. 

M = cor(genomes_SS)
M["w_max_2",] %>% sort()
M["D_6",] %>% sort()

D = as.data.frame(t(M["D_6",]))
D_cor = pivot_longer(D,H_1:h123_11)

D_cor = D_cor %>% mutate(SS_type = case_when(stringr::str_extract(name, "^.{2}") %in% c("H_", "D_") ~ 'SFS',
                           stringr::str_extract(name, "^.{2}") %in% c("h1", "h2") ~ 'Haplotype',
                           stringr::str_extract(name, "^.{3}") %in% c("Zns", "LD_","w_m") ~ 'LD',
                           TRUE ~ 'Other'))

ggplot(data = D_cor %>%
         dplyr::filter(SS_type == "LD"),
       aes(x = name, y = value)) +
  geom_point()
