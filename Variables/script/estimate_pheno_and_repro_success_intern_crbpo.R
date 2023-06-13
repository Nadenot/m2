
# -------------------------------------------------------------------------
# Title : estimate_pheno_and_repro_success_intern_crbpo
# Author : Paul Cuchot
# Date : 10/02/2023
# -------------------------------------------------------------------------

# the aim of this script is to estimate pheno and asymptot
# per year and per species 

# Library -----------------------------------------------------------------

# bayesian modeling + plot
library(R2jags)
library(MASS)
library(mcmcplots)

# general
library(tidyverse)
library(data.table)


# Load data ---------------------------------------------------------------

data <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")
# data <- read.csv("~/M2 EFCE/Stage/STOC_reporting-master/STOC_reporting-master/data_DB/data_reg.csv")



# restructure age 
# 0 = juvenile 
# 1 = adult
data <- data %>% mutate(AGE = ifelse(AGE == "1A", 0, 1))

## Trier les espèces de la + capturée à la - capturée et calculer les indicateurs phéno pour toutes les espèces jusqu'à ce que ça ne converge plus

# Tableau avec le nombre de captures par espèce  
sp_capt <- data %>% 
  group_by(ESPECE) %>% 
  count() %>% 
  arrange(by = desc(n))

# Add latin name to save file
sp_names <- fread("C:/git/STOC_reporting-master/library/nom_especes.csv")
sp_capt <- sp_capt %>% 
  left_join(sp_names[,1:2], by = c("ESPECE" = "SP"))
# write.csv(sp_capt, file = "C:/git/STOC/Variables/data/sp_capt.csv", row.names = FALSE)

species <- sp_capt$ESPECE

# # Test sur un jeu de données plus petit: SYLATR

data0 <- data
data <-  data <- data %>% 
  # filter(ESPECE == species[1])
  filter(ESPECE == "LOCLUS")

# function to summarize data per capture session
prod_sess <- function(data){
  prod <- data%>%
    group_by(YEAR, AGE, ID_PROG, JULIANDAY, ESPECE)%>%
    count() %>% # counts observations by group and put them in a column n
    pivot_wider(names_from = AGE, 
                values_from = n, 
                names_prefix = "n")%>%
    mutate(  
      n0 = ifelse(is.na(n0), 0, n0),
      n1 = ifelse(is.na(n1), 0, n1),
      prod = n0 / (n0+n1))%>%
    drop_na(prod)
  
  return(prod)
}


prod_f <- prod_sess(data)


# Structure data ----------------------------------------------------------


prod_f$an_sp <- paste(prod_f$YEAR, prod_f$ESPECE, sep = "_")

prod_f <- prod_f %>%
  mutate(an_f = as.factor(YEAR),
         an_n = as.numeric(an_f),
         site_f = as.factor(ID_PROG),
         site_n = as.numeric(site_f), # JAGS does not recognize years so we have to give them numbers
         sp_f = as.factor(ESPECE),
         sp_n = as.numeric(sp_f),  
         an_sp = as.factor(paste(an_f,sp_f, sep = "_")),
         an_sp_n = as.numeric(an_sp))

df_rd <- prod_f%>%
  distinct(an_sp_n, .keep_all = TRUE)

# data for the model
data <- list(nt = prod_f$n0+prod_f$n1,
             n0 = prod_f$n0,
             date = as.numeric(prod_f$JULIANDAY),
             N = nrow(prod_f),
             N_an = length(unique(prod_f$an_n)),
             N_site = length(unique(prod_f$site_n)),
             an = prod_f$an_n,
             sp = as.numeric(as.factor(as.character(prod_f$sp_n))),
             an2 = df_rd$an_n,
             sp2 = as.numeric(as.factor(as.character(df_rd$sp_n))),
             site = df_rd$site_n,
             N_sp = n_distinct(df_rd$ESPECE),
             N_an_sp_n = n_distinct(prod_f$an_sp_n),
             an_sp_n = prod_f$an_sp_n,
             rd_an_sp_n = df_rd$an_sp_n
)

# Model 1 -----------------------------------------------------------------

sink("model_crbpo")
cat("
    model{

  # loop on capture session

  for(i in 1:N){ # session data

    ## likelihood
    n0[i] ~ dbin(p[i], nt[i])
    
    p[i] <- asig[an_sp_n[i]]/(1+exp((csig[an_sp_n[i]]-date[i])/dsig[an_sp_n[i]]))
  
  }


  # loop on year and species

  for (ii in 1:N_an_sp_n){
    

  # csig parameter
  
    csig[ii] ~ dnorm(mu[ii], tau_res_csig[sp2[ii]]) 

    mu[ii] <- c[ii] + random_csig_site[site[ii]]  

    
  # scale parameter
  
    dsig[ii] ~ dnorm(mu_dsig[ii], tau_res_dsig[sp2[ii]])
    
    mu_dsig[ii] <- d[ii] + random_dsig_site[site[ii]]


  # asymptote parameter  
  
    asig[ii] ~ dnorm(mu_asig[ii], tau_res_asig[sp2[ii]])T(0.3,1)
    
    mu_asig[ii] <- a[ii] + random_asig_site[site[ii]]
    
     
  c[ii] ~ dnorm(150,0.01)

  a[ii] ~ dnorm(0,0.01)T(0,1)

  d[ii] ~ dnorm(0,0.01)T(0,10)
    
    
  }
  

  # random effect site
  
  for(z in 1:N_site){ #number of sites
 
  # random site 
    random_csig_site[z] ~ dnorm(0, tau_csig_site)
    random_asig_site[z] ~ dnorm(0, tau_asig_site)
    random_dsig_site[z] ~ dnorm(0, tau_dsig_site)
    
  }
  
  for(s in 1:N_sp){ # for  variance estimation (per species)
  
    sigma_res_csig[s] ~ dt(0, 0.01, 1)T(0,200) # Residual standard deviation
    sigma_res_dsig[s] ~ dt(0, 0.01, 1)T(0,10) # Residual standard deviation
    sigma_res_asig[s] ~ dt(0, 0.01, 1)T(0,1) # Residual standard deviation
  
    tau_res_csig[s] <- 1/(sigma_res_csig[s]*sigma_res_csig[s])
    tau_res_dsig[s] <- 1/(sigma_res_dsig[s]*sigma_res_dsig[s])
    tau_res_asig[s] <- 1/(sigma_res_asig[s]*sigma_res_asig[s])
  
  }
  
  # csig site
    sigma_csig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_csig_site <- pow(sigma_csig_site, -2)

  # dsig site
    sigma_dsig_site ~ dt(0, 0.01, 1)T(0,20)
    tau_dsig_site <- pow(sigma_dsig_site, -2)

  # asig site
    sigma_asig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_asig_site <- pow(sigma_asig_site, -2)
  
  

}", fill=TRUE)
sink()


# Initial values ----------------------------------------------------------


# init_f <-  function(){
#   list(a = rnorm(data$N_an_sp_n,150,15),
#        d = rnorm(data$N_an_sp_n,0,1),
#        c = rnorm(data$N_an_sp_n,0,1),
#        
#        random_csig_site = rep(0,data$N_site),
#        sigma_csig_site = runif(1,0,100),
#        
#        random_asig_site = rep(0,data$N_site),
#        sigma_asig_site = runif(1,0,1),
#        
#        random_dsig_site = rep(0,data$N_site),
#        sigma_dsig_site = runif(1,0,15)
#        
#   )
# }

# inits <- list(init1 = init_f(),init2=init_f(), init3=init_f())



# Run models --------------------------------------------------------------

# 1
parameters <- c("asig","csig","dsig",
                
                "d","a","c",
                
                "sigma_csig_site",
                "sigma_dsig_site",
                "sigma_asig_site",
                "random_csig_site",
                "random_dsig_site",
                "random_asig_site", "deviance"
)


md_1 <- jags(data = data,
             parameters.to.save = parameters,
             model.file = "model_crbpo",
             # inits = inits,
             n.chains = 3,
             n.iter = 10000,
             n.burnin = 3000)

## Il faut n.eff > 100 pour chaque paramètre (taille efficace, prend en compte l'autocorrélation)
## Rhat = Potential scale reduction factor = ratio entre la variabilité inter et intra chaîne 
# (si différence ça veut dire que les chaînes n'ont pas convergé)
# On considère comme acceptable Rhat < 1.1
sm <- md_1$BUGSoutput$summary
# only keep phenology estimates (c)
sm <- sm %>% 
  as.data.frame() %>% 
  mutate(parameter = rownames(sm)) %>% 
  dplyr::filter(str_detect(parameter, fixed("c["))) %>% 
  mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
  left_join(prod_f[,c("ID_PROG", "YEAR","ESPECE", "an_sp_n")])%>%
  distinct(mean,an_sp_n,ID_PROG,YEAR,ESPECE, .keep_all = TRUE) 


  

#saveRDS(md_1,"md_1.rds)


# library plot ------------------------------------------------------------

library(tidyverse)
library(bayesplot)

mat_md <- as.matrix(as.mcmc(md_1))


# jointure tableau pour identification espèce année
df_density <- mat_md%>%
  as.data.frame()%>%
  # select(starts_with("csig"))%>%
  dplyr::select(starts_with("c["))%>%
  pivot_longer(everything(),
               names_to = "parameter", values_to = "est")%>%
  mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
  left_join(prod_f[,c("ID_PROG", "YEAR","ESPECE", "an_sp_n")])%>%
  distinct(est,an_sp_n,ID_PROG,YEAR,ESPECE, .keep_all = TRUE) 



# plot des densités, un peu illisible
df_density%>%
  ggplot(aes(x = est,fill = parameter, group = parameter))+
  geom_density(alpha = 0.5)+
  facet_grid(ESPECE~YEAR, 
             switch = "y")+
  theme_minimal() +
  theme(
    strip.text.y.left = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  )


# get mean value for c (per species and per year)
df_density_mean <- df_density %>%
  group_by(YEAR, ESPECE)%>%
 # filter(YEAR < 2020)%>%
  summarize(mean_c = mean(est),
            sd__c = sd(est))

# evolution c par année (plot)
df_density%>%
  group_by(YEAR, ESPECE)%>%
  summarize(mean_c = mean(est),
            sd__c = sd(est))%>%
  ggplot(aes(x = YEAR, y = mean_c, color = ESPECE))+
  geom_point()+
  theme_bw()

# boxplot estimates per year and per species
df_density%>%
  ggplot(aes(x = as.factor(YEAR), y = est, color = ESPECE))+
  geom_boxplot()+
  xlab("Year") + ylab("xmid") +
  theme_bw()


################################################################################
#########################        SPECIES        ################################
################################################################################
data <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")
data <- data %>% mutate(AGE = ifelse(AGE == "1A", 0, 1))

### I will now implement a loop to create phenological indicators for all species as long as no n.eff <100
# With this condition, we stop after 4 species
# To be less restrictive, we can say that we want <10% of years that don't converge
# With this new criteria, we have 19 species
data0 <- data
conv <- TRUE #if models converges
i <- 1
while(conv){
  data_sp <- data0 %>% 
    dplyr::filter(ESPECE == species[i])
  prod_f <- prod_sess(data_sp)
    # Structure data ----------------------------------------------------------
  prod_f$an_sp <- paste(prod_f$YEAR, prod_f$ESPECE, sep = "_")
  prod_f <- prod_f %>%
  mutate(an_f = as.factor(YEAR),
           an_n = as.numeric(an_f),
           site_f = as.factor(ID_PROG),
           site_n = as.numeric(site_f), # JAGS does not recognize years so we have to give them numbers
           sp_f = as.factor(ESPECE),
           sp_n = as.numeric(sp_f),  
           an_sp = as.factor(paste(an_f,sp_f, sep = "_")),
           an_sp_n = as.numeric(an_sp))
  
  df_rd <- prod_f %>%
  distinct(an_sp_n, .keep_all = TRUE)
  
  # data for the model
  data <- list(nt = prod_f$n0+prod_f$n1,
               n0 = prod_f$n0,
               date = as.numeric(prod_f$JULIANDAY),
               N = nrow(prod_f),
               N_an = length(unique(prod_f$an_n)),
               N_site = length(unique(prod_f$site_n)),
               an = prod_f$an_n,
               sp = as.numeric(as.factor(as.character(prod_f$sp_n))),
               an2 = df_rd$an_n,
               sp2 = as.numeric(as.factor(as.character(df_rd$sp_n))),
               site = df_rd$site_n,
               N_sp = n_distinct(df_rd$ESPECE),
               N_an_sp_n = n_distinct(prod_f$an_sp_n),
               an_sp_n = prod_f$an_sp_n,
               rd_an_sp_n = df_rd$an_sp_n
  )
  
  # Model 1 -----------------------------------------------------------------
  
  sink("model_crbpo")
  cat("
    model{
  # loop on capture session
  for(i in 1:N){ # session data
    ## likelihood
    n0[i] ~ dbin(p[i], nt[i])
    p[i] <- asig[an_sp_n[i]]/(1+exp((csig[an_sp_n[i]]-date[i])/dsig[an_sp_n[i]]))
  }
  # loop on year and species
  for (ii in 1:N_an_sp_n){
  # csig parameter
    csig[ii] ~ dnorm(mu[ii], tau_res_csig[sp2[ii]]) 
    mu[ii] <- c[ii] + random_csig_site[site[ii]]  
  # scale parameter
    dsig[ii] ~ dnorm(mu_dsig[ii], tau_res_dsig[sp2[ii]])
    mu_dsig[ii] <- d[ii] + random_dsig_site[site[ii]]
  # asymptote parameter  
    asig[ii] ~ dnorm(mu_asig[ii], tau_res_asig[sp2[ii]])T(0.3,1)
    mu_asig[ii] <- a[ii] + random_asig_site[site[ii]]
  c[ii] ~ dnorm(150,0.01)
  a[ii] ~ dnorm(0,0.01)T(0,1)
  d[ii] ~ dnorm(0,0.01)T(0,10)
  }
  # random effect site
  for(z in 1:N_site){ #number of sites
  # random site 
    random_csig_site[z] ~ dnorm(0, tau_csig_site)
    random_asig_site[z] ~ dnorm(0, tau_asig_site)
    random_dsig_site[z] ~ dnorm(0, tau_dsig_site)
  }
  for(s in 1:N_sp){ # for  variance estimation (per species)
    sigma_res_csig[s] ~ dt(0, 0.01, 1)T(0,200) # Residual standard deviation
    sigma_res_dsig[s] ~ dt(0, 0.01, 1)T(0,10) # Residual standard deviation
    sigma_res_asig[s] ~ dt(0, 0.01, 1)T(0,1) # Residual standard deviation
    tau_res_csig[s] <- 1/(sigma_res_csig[s]*sigma_res_csig[s])
    tau_res_dsig[s] <- 1/(sigma_res_dsig[s]*sigma_res_dsig[s])
    tau_res_asig[s] <- 1/(sigma_res_asig[s]*sigma_res_asig[s])
  }
  # csig site
    sigma_csig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_csig_site <- pow(sigma_csig_site, -2)
  # dsig site
    sigma_dsig_site ~ dt(0, 0.01, 1)T(0,20)
    tau_dsig_site <- pow(sigma_dsig_site, -2)
  # asig site
    sigma_asig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_asig_site <- pow(sigma_asig_site, -2)
}", fill=TRUE)
  sink()
   # Run models --------------------------------------------------------------
  parameters <- c("asig","csig","dsig",
                  
                  "d","a","c",
                  
                  "sigma_csig_site",
                  "sigma_dsig_site",
                  "sigma_asig_site",
                  "random_csig_site",
                  "random_dsig_site",
                  "random_asig_site", "deviance"
  )
  
  md_1 <- jags(data = data,
               parameters.to.save = parameters,
               model.file = "model_crbpo",
               # inits = inits,
               n.chains = 3,
               n.iter = 70000,
               n.burnin = 20000)
  
  ## Il faut n.eff > 100 pour chaque paramètre (taille efficace, prend en compte l'autocorrélation)
  ## Rhat = Potential scale reduction factor = ratio entre la variabilité inter et intra chaîne 
  # (si différence ça veut dire que les chaînes n'ont pas convergé)
  # On considère comme acceptable Rhat < 1.1
  sm <- md_1$BUGSoutput$summary
  # only keep phenology estimates (c)
  sm <- sm %>% 
    as.data.frame() %>% 
    mutate(names = rownames(sm)) %>% 
    dplyr::filter(str_detect(names, fixed("c[")))
  neff100 <- sm %>% dplyr::filter(n.eff < 100)
  rhat11 <- sm %>% dplyr::filter(Rhat > 1.1)  
  if(nrow(neff100)>0 | nrow(rhat11)>0){ # if more than 3 years (among 35 so 10%) don't converge, we stop
    conv <- FALSE
  } else i <- i+1 #else we continue with the next species
  
  # Save estimates
  mat_md <- as.matrix(as.mcmc(md_1))
    # jointure tableau pour identification espèce année
  df_density <- mat_md%>%
    as.data.frame()%>%
    # select(starts_with("csig"))%>%
    dplyr::select(starts_with("c["))%>%
    pivot_longer(everything(),
                 names_to = "parameter", values_to = "est")%>%
    mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
    left_join(prod_f[,c("ID_PROG", "YEAR","ESPECE", "an_sp_n")])%>%
    distinct(est,an_sp_n,ID_PROG,YEAR,ESPECE, .keep_all = TRUE)
  
  # get mean value for c (per species and per year)
  df_density_mean <- df_density %>%
    group_by(YEAR, ESPECE)%>%
    summarize(mean_c = mean(est),
              sd__c = sd(est))
  
  if(i == 2){ #first round but we already did i <- i+1
    df_pheno <- df_density
    df_pheno_mean <- df_density_mean
  } 
  if(i>2){
    df_pheno <- rbind(df_pheno, df_density)
    df_pheno_mean <- rbind(df_pheno_mean, df_density_mean)
  }
}


write.csv(df_pheno_mean, file = "C:/git/STOC/Variables/data/pheno.csv", row.names = FALSE)
write.csv(df_pheno, file = "C:/git/STOC/Variables/data/pheno_detail.csv", row.names = FALSE)


# graphs to check convergence
traceplot(md_1, varname = "c")



################################################################################
########################        50 SPECIES        ##############################
################################################################################

# Here I don't want to stop when it stops converging, I want to test the quality of the output values
# So I decide to consider 50 species and test Rhat and n.eff afterwards

data <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")
data <- data %>% mutate(AGE = ifelse(AGE == "1A", 0, 1))

## Trier les espèces de la + capturée à la - capturée et calculer les indicateurs phéno pour toutes les espèces jusqu'à ce que ça ne converge plus

sp_capt <- fread("C:/git/STOC/Variables/data/sp_capt.csv")

species <- sp_capt$ESPECE

data0 <- data

for(i in 1:50){
  data_sp <- data0 %>% 
    dplyr::filter(ESPECE == species[i])
  prod_f <- prod_sess(data_sp)
  # Structure data ----------------------------------------------------------
  prod_f$an_sp <- paste(prod_f$YEAR, prod_f$ESPECE, sep = "_")
  prod_f <- prod_f %>%
    mutate(an_f = as.factor(YEAR),
           an_n = as.numeric(an_f),
           site_f = as.factor(ID_PROG),
           site_n = as.numeric(site_f), # JAGS does not recognize years so we have to give them numbers
           sp_f = as.factor(ESPECE),
           sp_n = as.numeric(sp_f),  
           an_sp = as.factor(paste(an_f,sp_f, sep = "_")),
           an_sp_n = as.numeric(an_sp))
  
  df_rd <- prod_f %>%
    distinct(an_sp_n, .keep_all = TRUE)
  
  # data for the model
  data <- list(nt = prod_f$n0+prod_f$n1,
               n0 = prod_f$n0,
               date = as.numeric(prod_f$JULIANDAY),
               N = nrow(prod_f),
               N_an = length(unique(prod_f$an_n)),
               N_site = length(unique(prod_f$site_n)),
               an = prod_f$an_n,
               sp = as.numeric(as.factor(as.character(prod_f$sp_n))),
               an2 = df_rd$an_n,
               sp2 = as.numeric(as.factor(as.character(df_rd$sp_n))),
               site = df_rd$site_n,
               N_sp = n_distinct(df_rd$ESPECE),
               N_an_sp_n = n_distinct(prod_f$an_sp_n),
               an_sp_n = prod_f$an_sp_n,
               rd_an_sp_n = df_rd$an_sp_n
  )
  
  # Model 1 -----------------------------------------------------------------
  
  sink("model_crbpo")
  cat("
    model{
  # loop on capture session
  for(i in 1:N){ # session data
    ## likelihood
    n0[i] ~ dbin(p[i], nt[i])
    p[i] <- asig[an_sp_n[i]]/(1+exp((csig[an_sp_n[i]]-date[i])/dsig[an_sp_n[i]]))
  }
  # loop on year and species
  for (ii in 1:N_an_sp_n){
  # csig parameter
    csig[ii] ~ dnorm(mu[ii], tau_res_csig[sp2[ii]]) 
    mu[ii] <- c[ii] + random_csig_site[site[ii]]  
  # scale parameter
    dsig[ii] ~ dnorm(mu_dsig[ii], tau_res_dsig[sp2[ii]])
    mu_dsig[ii] <- d[ii] + random_dsig_site[site[ii]]
  # asymptote parameter  
    asig[ii] ~ dnorm(mu_asig[ii], tau_res_asig[sp2[ii]])T(0.3,1)
    mu_asig[ii] <- a[ii] + random_asig_site[site[ii]]
  c[ii] ~ dnorm(150,0.01)
  a[ii] ~ dnorm(0,0.01)T(0,1)
  d[ii] ~ dnorm(0,0.01)T(0,10)
  }
  # random effect site
  for(z in 1:N_site){ #number of sites
  # random site 
    random_csig_site[z] ~ dnorm(0, tau_csig_site)
    random_asig_site[z] ~ dnorm(0, tau_asig_site)
    random_dsig_site[z] ~ dnorm(0, tau_dsig_site)
  }
  for(s in 1:N_sp){ # for  variance estimation (per species)
    sigma_res_csig[s] ~ dt(0, 0.01, 1)T(0,200) # Residual standard deviation
    sigma_res_dsig[s] ~ dt(0, 0.01, 1)T(0,10) # Residual standard deviation
    sigma_res_asig[s] ~ dt(0, 0.01, 1)T(0,1) # Residual standard deviation
    tau_res_csig[s] <- 1/(sigma_res_csig[s]*sigma_res_csig[s])
    tau_res_dsig[s] <- 1/(sigma_res_dsig[s]*sigma_res_dsig[s])
    tau_res_asig[s] <- 1/(sigma_res_asig[s]*sigma_res_asig[s])
  }
  # csig site
    sigma_csig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_csig_site <- pow(sigma_csig_site, -2)
  # dsig site
    sigma_dsig_site ~ dt(0, 0.01, 1)T(0,20)
    tau_dsig_site <- pow(sigma_dsig_site, -2)
  # asig site
    sigma_asig_site ~ dt(0, 0.01, 1)T(0,1)
    tau_asig_site <- pow(sigma_asig_site, -2)
}", fill=TRUE)
  sink()
  # Run models --------------------------------------------------------------
  parameters <- c("asig","csig","dsig",
                  "d","a","c",
                  "sigma_csig_site",
                  "sigma_dsig_site",
                  "sigma_asig_site",
                  "random_csig_site",
                  "random_dsig_site",
                  "random_asig_site", "deviance"
  )
  
  md_1 <- jags(data = data,
               parameters.to.save = parameters,
               model.file = "model_crbpo",
               # inits = inits,
               n.chains = 3,
               n.iter = 70000,
               n.burnin = 20000)
  
  
  ## Il faut n.eff > 100 pour chaque paramètre (taille efficace, prend en compte l'autocorrélation)
  ## Rhat = Potential scale reduction factor = ratio entre la variabilité inter et intra chaîne 
  # (si différence ça veut dire que les chaînes n'ont pas convergé)
  # On considère comme acceptable Rhat < 1.1
  sm <- md_1$BUGSoutput$summary
  # only keep phenology estimates (c)
  sm <- sm %>% 
    as.data.frame() %>% 
    mutate(parameter = rownames(sm)) %>% 
    dplyr::filter(str_detect(parameter, fixed("c["))) %>% 
    mutate(an_sp_n = as.numeric(gsub(".*?([0-9]+).*", "\\1", parameter)))%>%
    left_join(prod_f[,c("YEAR","ESPECE", "an_sp_n")])%>%
    distinct(mean,an_sp_n,YEAR,ESPECE, .keep_all = TRUE) 
  
  if(i == 1){ 
    df_pheno_conv <- sm
  } 
  if(i>1){
    df_pheno_conv <- rbind(df_pheno_conv, sm)
  }
}


write.csv(df_pheno_conv, file = "C:/git/STOC/Variables/data/pheno_conv_model.csv", row.names = FALSE)







################################################################################
################      Plot phenology for each species      #####################
################################################################################

library(ggplot2)
# pheno <- fread("C:/git/STOC/Variables/data/pheno.csv")
pheno <- fread("C:/git/STOC/Variables/data/pheno33.csv")
df_pheno_conv <- fread("C:/git/STOC/Variables/data/pheno_conv_model.csv")

# sort by pheno

pheno_sorted <- df_pheno_conv %>% arrange(by = mean) %>% mutate(ESPECE = as.factor(ESPECE)) %>% 
  mutate(Date = as.Date(round(mean), origin = as.Date("2023-01-01"))) %>% 
  rename(mean_c = mean)

pheno_sorted <- pheno %>% arrange(by = mean_c) %>% mutate(ESPECE = as.factor(ESPECE)) %>% 
  mutate(Date = as.Date(round(mean_c), origin = as.Date("2023-01-01"))) 

gg <- ggplot(pheno_sorted, aes(x=reorder(ESPECE, mean_c), y=Date)) + 
  geom_boxplot()
gg <- gg + scale_y_date(labels = scales::date_format("%d/%m"))
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + ggtitle("Variabilité interannuelle de la date d'envol moyenne") +
  xlab("Espèce") + ylab("Date d'envol moyenne")
gg

# Group species by migration status  
# Add migration  
data_mod <- fread("C:/git/STOC/Variables/data/model/data_model.csv")
pheno_sorted <- pheno_sorted %>% left_join(data_mod[,c("ESPECE", "MIGRATION")]) %>% distinct()
# plot
gg <- ggplot(pheno_sorted, aes(x=reorder(ESPECE, mean_c), y=Date, fill = MIGRATION)) + 
  geom_boxplot()
gg <- gg + scale_y_date(labels = scales::date_format("%d/%m"))
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + ggtitle("Variabilité interannuelle de la date d'envol moyenne") +
  xlab("Espèce") + ylab("Date d'envol moyenne")
gg

# French species name  
names <- fread("C:/git/STOC_reporting-master/library/nom_especes.csv", encoding = "Latin-1") 
names <- names %>% rename(ESPECE = SP) %>% dplyr::select(ESPECE, nom_sc, nom_fr)
pheno_sorted <- pheno_sorted %>% left_join(names)
# plot  
gg <- ggplot(pheno_sorted, aes(x=reorder(nom_fr, mean_c), y=Date, fill = MIGRATION)) + 
  geom_boxplot()
gg <- gg + scale_y_date(labels = scales::date_format("%d/%m"))
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + ggtitle("Variabilité interannuelle de la date d'envol moyenne") +
  xlab("Espèce") + ylab("Date d'envol moyenne")
gg


################################################################################
###############   Control quality of phenology indicator   #####################
################################################################################


################     Number of sites capturing species     #####################

data <- fread("C:/git/STOC_reporting-master/data_DB/data.csv")
pheno <- pheno_sorted
# Count number of captured individuals per species  
capt_ind_site <-data %>% 
  filter(ESPECE %in% unique(pheno$ESPECE)) %>% 
  group_by(YEAR, ESPECE) %>% 
  count() %>% 
  rename(tot_capt_ind = n)

# Count mean number of captured individuals per site  
capt_ind_site_m <-data %>% 
  filter(ESPECE %in% unique(pheno$ESPECE)) %>% 
  group_by(NEW.ID_PROG, YEAR, ESPECE) %>% 
  count() %>% 
  group_by(YEAR, ESPECE) %>% 
  summarise(mean_capt_site = round(mean(n, na.rm = TRUE)))

# Count sites who captured each species
nb_site_sp <-data %>% 
  filter(ESPECE %in% unique(pheno$ESPECE)) %>% 
  group_by(NEW.ID_PROG, YEAR, ESPECE) %>% 
  count() %>% 
  group_by(ESPECE, YEAR) %>% 
  count() %>% 
  rename(nb_sites = n)

# merge with pheno
data_plot <- pheno_sorted %>% 
  left_join(capt_ind_site) %>% 
  left_join(capt_ind_site_m) %>% 
  left_join(nb_site_sp)


### plot ###
# Total number of captured individuals 
sp_capt <- fread("C:/git/STOC/cluster/data/sp_capt.csv")
species <- sp_capt[1:50,1]

# delete within session recaptures
data <- data %>% 
  group_by(DATE, BAGUE) %>% 
  slice(1) %>% # takes the first occurrence if there is a tie
  ungroup()

capt_ind_site <- data %>% 
  filter(ESPECE %in% species$ESPECE) %>%
  filter(!(ESPECE %in% c("PANBIA", "SYLALA"))) %>% 
  filter(YEAR > 2000) %>% 
  group_by(YEAR, ESPECE) %>% 
  count() %>% 
  rename(tot_capt_ind = n)

# French species name  
names <- fread("C:/git/STOC_reporting-master/library/nom_especes.csv", encoding = "Latin-1") 
names <- names %>% rename(ESPECE = SP) %>% dplyr::select(ESPECE, nom_sc, nom_fr)
capt_ind_site <- capt_ind_site %>% left_join(names)

gg <- ggplot(capt_ind_site, aes(x=reorder(nom_sc, tot_capt_ind), y=tot_capt_ind)) + 
  geom_boxplot()
gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Number of individuals captured each year", subtitle = "Interannual variability") +
  xlab("Species") + ylab("Number of captured individuals")
gg


# plot nb sites
gg <- ggplot(data_plot, aes(x=reorder(nom_fr, mean_c), y=nb_sites)) + 
  geom_boxplot()
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Nombre de stations capturant chaque espèce", subtitle = "Variabilité interannuelle") +
  xlab("Espèce") + ylab("Nombre de stations")
gg

# plot mean captured individuals / site
gg <- ggplot(data_plot, aes(x=reorder(nom_fr, mean_c), y=mean_capt_site)) + 
  geom_boxplot()
gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Nombre moyen d'individus capturés par site par espèce", subtitle = "Variabilité interannuelle") +
  xlab("Espèce") + ylab("Nombre d'individus capturés par an (moyenne des sites)")
gg

# there were very few stations at the beginning so we will see what it looks like of we consider only years > 2005  
data_plot2005 <-  data_plot %>% filter(YEAR > 2005)

## plot nb sites
gg <- ggplot(data_plot2005, aes(x=reorder(nom_fr, nb_sites), y=nb_sites)) + 
  geom_boxplot()
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Nombre de stations capturant chaque espèce", subtitle = "Variabilité interannuelle depuis 2005") +
  xlab("Espèce") + ylab("Nombre de stations")
gg

## plot mean captured individuals / site
gg <- ggplot(data_plot, aes(x=reorder(nom_fr, mean_c), y=mean_capt_site)) + 
  geom_boxplot()
gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Nombre moyen d'individus capturés par site par espèce", subtitle = "Variabilité interannuelle") +
  xlab("Espèce") + ylab("Nombre d'individus capturés par an (moyenne des sites)")
gg


###################        Get confidence interval           ###################

# 95% confidence interval from sd and sample size
pheno <- fread("C:/git/STOC/Variables/data/pheno33.csv")

# Compute the size
pheno_int <- data_plot %>% 
  mutate(se = sd__c / sqrt(nb_sites)) %>% #se = standard error ; n = number of sites capturing each species each year
  mutate(low_confint = mean_c - 1.96 * se) %>% 
  mutate(up_confint = mean_c + 1.96 * se) %>% 
  mutate(int = round(up_confint - low_confint, 1))

# plot confidence interval 
pheno_int2000 <- pheno_int %>% filter(YEAR > 1999)
# plot
gg <- ggplot(pheno_int2000, aes(x=reorder(nom_fr, mean_c), y=int)) + 
  geom_boxplot()
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "95% confidence interval", subtitle = "Variabilité interannuelle depuis 2000") +
  xlab("Espèce") + ylab("largeur de l'intervalle (jours)")
gg


# 95% quantile  
pheno_detail <- fread("C:/git/STOC/Variables/data/pheno_detail.csv")
# test with PARMAJ 2014
# with the 95% conf int, we get 136.2959 - 137.2316 --> 0.9
# test <- pheno_detail %>% filter(ESPECE == "PARMAJ", YEAR == 2014)
# quantile(test$est, probs = 0.025, na.rm = TRUE) # 132.2487  
# quantile(test$est, probs = 0.975, na.rm = TRUE) # 142.1563  
pheno_quantile <- pheno_detail %>% 
  group_by(YEAR, ESPECE)%>%
  summarize(mean_c = mean(est),
            quant025 = quantile(est, probs = 0.025, na.rm = TRUE),
            quant975 = quantile(est, probs = 0.975, na.rm = TRUE)) %>% 
  mutate(int = quant975 - quant025)
# plot
names <- fread("C:/git/STOC_reporting-master/library/nom_especes.csv", encoding = "Latin-1") 
names <- names %>% rename(ESPECE = SP) %>% dplyr::select(ESPECE, nom_sc, nom_fr)
pheno_quantile <- pheno_quantile %>% left_join(names)
pheno_quantile2000 <- pheno_quantile %>% filter(YEAR > 1999)

gg <- ggplot(pheno_quantile2000, aes(x=reorder(nom_fr, mean_c), y=int)) + 
  geom_boxplot()
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "95% quantile (95% des données dans cet intervalle)", subtitle = "Variabilité interannuelle depuis 2000") +
  xlab("Espèce") + ylab("largeur de l'intervalle (jours)")
gg

data_plot_2000 <- data_plot_ %>% filter(YEAR >1999)
gg <- ggplot(data_plot_2000, aes(x=reorder(nom_fr, tot_capt_ind), y=int95)) + 
  geom_boxplot()
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "95% quantile (95% des données dans cet intervalle)", subtitle = "Variabilité interannuelle depuis 2000") +
  xlab("Espèce (de la plus à la moins capturée)") + ylab("largeur de l'intervalle (jours)")
gg



#### Get % of years without convergence
sm <- fread("C:/git/STOC/Variables/data/pheno_conv_model.csv")

sm_conv_years <- sm %>% 
  mutate(conv = ifelse(Rhat < 1.1 & n.eff > 100, TRUE, FALSE)) %>% 
  group_by(ESPECE, conv) %>% 
  count() %>% 
  pivot_wider(names_from = conv, values_from = n) %>% 
  


#######   Relationship between pheno estimation quality and sample size  #######
data_plot_ <- data_plot %>% 
  left_join(sm) %>% 
  dplyr::filter(YEAR>1999)

data_mod <- fread("C:/git/STOC/Variables/data/model/data_model_all_a.csv")
data_mod_pheno <- data_mod %>% 
  dplyr::select(ID_PROG:MIGRATION) %>% 
  left_join(data_plot_) %>% 
  filter(!is.na(Rhat))

## Rhat  
# ~ nb ind capt
gg <- ggplot(data_plot_, aes(x=tot_capt_ind, y=Rhat)) + geom_point()
gg
GAMM <- gamm(Rhat ~ s(tot_capt_ind), random = list(YEAR = ~1, ID_PROG = ~1), 
             family = poisson, data = data_mod_pheno)  
plot(GAMM$gam,pages=1)

 # ~ nb ind/site  
gg <- ggplot(data_plot_, aes(x=mean_capt_site, y=Rhat)) + geom_point()
gg

# ~ nb sites
gg <- ggplot(data_plot_, aes(x=nb_sites, y=Rhat)) + geom_point()
gg

## n.eff  
# ~ nb ind capt
gg <- ggplot(data_plot_, aes(x=tot_capt_ind, y=n.eff)) + geom_point()
gg
# ~ nb ind/site  
gg <- ggplot(data_plot_, aes(x=mean_capt_site, y=n.eff)) + geom_point()
gg
# ~ nb sites
gg <- ggplot(data_plot_, aes(x=nb_sites, y=n.eff)) + geom_point()
gg

## Variance  
## Sd
# ~ nb ind capt
gg <- ggplot(data_plot_, aes(x=tot_capt_ind, y=sd, colour = ESPECE)) + geom_point()
gg
# ~ nb ind/site  
gg <- ggplot(data_plot_, aes(x=mean_capt_site, y=sd, colour = ESPECE)) + geom_point()
gg
# ~ nb sites
gg <- ggplot(data_plot_, aes(x=nb_sites, y=sd, colour = ESPECE)) + geom_point()
gg

## 95% quantile
data_plot_ <- data_plot_ %>% 
  mutate(int95 = `97.5%` - `2.5%`)
# ~ nb ind capt
gg <- ggplot(data_plot_, aes(x=tot_capt_ind, y=int95, colour = ESPECE)) + geom_point() #+ theme(legend.position='none')
gg
# ~ nb ind/site  
gg <- ggplot(data_plot_, aes(x=mean_capt_site, y=int95, colour = ESPECE)) + geom_point()
gg
# ~ nb sites
gg <- ggplot(data_plot_, aes(x=nb_sites, y=int95, colour = ESPECE)) + geom_point()
gg

#### Variance of estimation precision  
data_plot_sd <- data_plot %>% 
  mutate(int95 = `97.5%` - `2.5%`)
# We have data for each year-species, and we want the sd of 95% interval  between years for each species  
data_plot_sd <- data_plot_sd %>% 
  group_by(nom_fr) %>% 
  summarise(sd_int = sd(int95), mean_int = mean(int95), mean_ind = mean(tot_capt_ind))

# plot
# hist  
hist(data_plot_sd$sd_int, breaks = 15, main = "Répartition écart-types de la variation \n interannuelle de l'intervalle 95%",
     xlab = "sd interannuel 95% interval")

gg <- ggplot(data_plot_sd, aes(x=reorder(nom_fr, mean_ind), y=sd_int)) + 
  geom_bar(stat = "identity")
# gg <- gg + scale_y_continuous(trans='log2')
gg <- gg + coord_flip() # Tourner le box plot
# gg <- gg + geom_jitter(position=position_jitter(0), colour = "#737373") # add points
gg <- gg + labs(title = "Ecart-type de la variation interannuelle de l'intervalle 95%") +
  xlab("Espèce (de la plus à la moins capturée)") + ylab("sd")
gg

gg <- ggplot(data_plot_sd, aes(x=mean_ind, y=sd_int)) + geom_point() #+ theme(legend.position='none')
# gg <- gg + scale_x_continuous(trans='log2')
gg <- gg +  xlab("sd interannuel 95% interval") + ylab("Nombre moyen d'indivdus capturés")
gg



## Correlation Moussus-Cuchot  
corr <- fread("C:/git/STOC/Variables/data/corr_decalage_pheno_moussus_cuchot.csv")
summary(corr)
sd(corr$correlation_MC)
