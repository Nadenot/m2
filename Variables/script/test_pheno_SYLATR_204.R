## Here I will see if I can see an evolution in the age of juveniles captured at PY's station in SYLATR
#I want to know if we first capture a lot of PUL, and then 1A, which would mean that they did not just fledge but may have dispersed

data <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")

# Select PY's station  
dataPY <- data %>% 
  filter(NEW.ID_PROG == 204) %>% 
  filter(ESPECE == "SYLATR") # focus on SYLATR

# when did they start using PUL (poussins) ?
min(dataPY$YEAR[which(dataPY$AGE == "PUL")])
# 2009

# Plot % of PUL at each session
# % total captures 
propPUL <- dataPY %>% 
  # keep categories PUL, 1A & group the rest as AD
  mutate(AGE2 = ifelse(AGE %in% c("PUL", "1A"),AGE, "AD" )) %>% 
  group_by(YEAR, DATE, AGE2) %>% 
  count() %>% 
  pivot_wider(names_from = AGE2, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(propPULtot = round(PUL/(PUL + `1A` + AD),1)) %>% 
  mutate(propPULjuv = ifelse(`1A` != 0, round(PUL/(PUL + `1A`),1), 0))

# count session
propPUL$session <- 1
for(i in 2:nrow(propPUL)){
  if(propPUL$YEAR[i] == propPUL$YEAR[i-1]) propPUL$session[i] <- propPUL$session[i-1] + 1
}


# Il y a que 5 années au total avec des captures notées PUL, et ce ne sont pas les années les plus récentes, 
# donc je ne pense pas que ce soit exploitable
gg <- ggplot(data = propPUL, aes(x = as.factor(session), y = propPULtot)) +
  geom_bar(stat="identity")
gg



################################################################################
data3session <- fread("C:/git/STOC_reporting-master/data_DB/data3session.csv")

# par contre je vais garde l'idée de plot le % de jeunes par session
data3session <- data3session %>% mutate(AGE = ifelse(AGE == "1A", "JUV", "AD"))

propJUV <- data3session %>% 
  group_by(ESPECE, YEAR, ID_PROG,DATE, AGE_first) %>% 
  count() %>% 
  pivot_wider(names_from = AGE_first, values_from = n) %>% 
  replace(is.na(.), 0) %>% 
  mutate(prod = round(JUV/(JUV + AD),1))


prodSYLATR <- propJUV %>% filter(ESPECE == "SYLATR")
# count session
prodSYLATR$session <- 1
for(i in 2:nrow(prodSYLATR)){
  if(prodSYLATR$YEAR[i] == prodSYLATR$YEAR[i-1] & prodSYLATR$ID_PROG[i] == prodSYLATR$ID_PROG[i-1]) prodSYLATR$session[i] <- prodSYLATR$session[i-1] + 1
}

# mean prod per session  
prodSYLATR_ <- prodSYLATR %>% 
  group_by(YEAR, session) %>% 
  summarise(mean_prod = mean(prod, na.rm = TRUE))
