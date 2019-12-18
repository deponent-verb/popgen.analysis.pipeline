#debugging NA values script

pacman::p_load(tidyverse)

df<-read_csv("./data/toy_df.csv")

#extract H values
temp<-df[,1:14]
temp2<-temp %>% filter_all(any_vars(is.na(.))) 


#extract D values
temp<-df[,15:24]
temp1<-temp %>% filter_all(any_vars(is.na(.))) 

#D is not outputting NAs. 

#new_data <- data %>% filter_all(any_vars(is.na(.))) 
