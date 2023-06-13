# Steps data selection  
library(data.table)

# raw data
dataSTOC <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")
# 460138 obs

# delete within session recatures
dataSTOC_recatch <- dataSTOC %>% 
  group_by(DATE, BAGUE) %>% 
  #filter(DATE == min(DATE)) %>% 
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()
# 414958 obs

