library(DataExplorer)

var_mod <- fread("C:/git/STOC/Variables/data/model/var_model_final_allsp_scaled.csv")

#normality 
plot_histogram(var_mod) 
plot_qq(var_mod) 

var_mod_num <- var_mod %>% select_if(is.numeric) %>% select(6:26) %>% 
  filter(!is.na(meanNDVI)) %>%
  filter(!is.na(early_SPEI)) %>% 
  filter(!is.na(early_nb_ECE_SPEI))

library(corrplot)

## All variables
M <- cor(var_mod_num)
corrplot(M, method = "number", type="upper")

## Variables in the best model  
var_mod_best <- var_mod %>% select_if(is.numeric) %>% select(7:14) %>% 
  filter(!is.na(meanNDVI)) %>%
  filter(!is.na(early_SPEI)) 

M <- cor(var_mod_best)
corrplot(M, method = "number", type="upper")

