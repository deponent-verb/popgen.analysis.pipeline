#example script for parallel coordinates plot

pacman::p_load(tidyverse,popgen.tools)

#load data
DF<-readRDS("~/work/MPhil/data/toy_data.rds")

#sanity checking the dataframe

#parallel coordinations plot

DF %>% 
  sample_n(2000) %>% 
  gather(key = "SS", value = "value", D:h123) %>% 
  filter(SS != "H") %>% 
  ggplot() + 
  aes(SS, value) + 
  geom_line(aes(group = ID), alpha = 0.1) + 
  geom_smooth(aes(group = sweep))

DF %>% 
  sample_n(2000) %>% 
  gather(key = "SS", value = "value", D:h123) %>% 
  filter(SS != "H") %>% 
  ggplot() + 
  aes(sweep, value,fill=position) + 
  geom_boxplot() + 
  facet_wrap(~SS, scales = "free") 
  


DF %>% 
  select(D:h123) %>% 
  map(~lm(.x ~ sweep + position, data = DF)) %>% 
  map_df(broom::tidy, .id = "Pred") %>% 
  filter(p.value <= 0.05) %>% 
  print(n = 30)



