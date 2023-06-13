require(tidyverse)
require(dplyr)
require(data.table)
require(ggplot2)
require(glmmTMB)
require(brms)
require(phytools)


## Full model per species  

cat(" ------------------    Import data    -----------------------")

var_model_scaled <- fread("/scratchbeta/adenotn/data/var_model_final_allsp_scaled.csv")
# var_model_scaled <- fread("C:/Users/Nathalie Adenot/STOC/cluster/data/var_model_final_allsp_scaled.csv")
# var_model_scaled <- fread("C:/git/STOC/cluster/data/var_model_final_allsp_scaled.csv")

# there are 162 species, we will run the model on the 50 most captured  

sp_capt <- fread("/scratchbeta/adenotn/data/sp_capt.csv")
# sp_capt <- fread("C:/Users/Nathalie Adenot/STOC/cluster/data/sp_capt.csv")
# sp_capt <- fread("C:/git/STOC/cluster/data/sp_capt.csv")
sp_capt <- sp_capt %>% mutate(phylo = gsub(" ", "_", nom_sc)) %>% 
  mutate(phylo = gsub("Cyanistes", "Parus", phylo)) # in the tree, Cyanistes caeruleus is Parus caeruleus



species <- sp_capt$ESPECE

# keep the 50 most captured species, without 
var_model_scaled50 <- var_model_scaled %>% filter(ESPECE %in% species[1:50]) %>% 
  filter(!(ESPECE %in% c("PANBIA","SYLALA")))



cat(" ------------------   RUN glmmTMB   -----------------------")

mod <- glmmTMB(cbind(JUV, AD) ~ mean_BAL + cat + meanNDVI + density +
                 early_SPEI + late_SPEI +
                 I(early_SPEI^2) + I(late_SPEI^2) +
                 mean_aT_early + mean_aT_late + I(mean_aT_early^2) + I(mean_aT_late^2) +
                 mean_aNDVI_early + mean_aNDVI_late + mean_aNDVI_early:meanNDVI +
                early_SPEI:ESPECE + late_SPEI:ESPECE + cat:early_SPEI + cat:late_SPEI +
                mean_BAL:I(early_SPEI^2) + mean_BAL:I(late_SPEI^2) +
                mean_aT_early:ESPECE + mean_aT_late:ESPECE + 
                cat:mean_aT_early + cat:mean_aT_late + I(mean_aT_early^2):cat + I(mean_aT_late^2):cat +
                 mean_aNDVI_early:ESPECE + mean_aNDVI_late:ESPECE +
                ESPECE + (1|ID_PROG) + (1|YEAR),
              family=binomial, data=var_model_scaled50)


s <- summary(mod)
coef <- as.data.frame(s[["coefficients"]][["cond"]])
coef$var <- row.names(coef)
# extract species name
coef_sp <- coef %>% filter(grepl("ESPECE", var) & nchar(var)>12) %>%
  mutate(Species = substr(var,nchar(var)-6+1, nchar(var)))
# extract variable name 
coef_sp <- coef_sp %>% 
  mutate(var_name = substr(var,1, nchar(var)-13))

# mod <- glmmTMB(cbind(JUV, AD) ~ mean_BAL + cat + meanNDVI + density +
#                 early_nb_ECE_aTemp + late_nb_ECE_aTemp + cat:early_nb_ECE_aTemp + cat:late_nb_ECE_aTemp +
#                 early_SPEI + late_nb_ECE_SPEI +
#                 mean_aNDVI_early + meanNDVI:mean_aNDVI_early +
#                 (1|ID_PROG) + (1|YEAR),
#               family=binomial, data=var_model_scaled_noNA_ECE)





cat(" ------------------   Save file   -----------------------")

write.csv(coef_sp, file = "/scratchbeta/adenotn/output/coef_sp_meta2.csv", row.names = FALSE)
# write.csv(coef_sp, file = "C:/git/STOC/stats/coef_sp_meta2.csv", row.names = FALSE)




cat(" ------------------    Prepare data brms    -----------------------")

cat(" ------------------    Import tree    -----------------------")
phylo <- ape::read.nexus("/scratchbeta/adenotn/data/AllBirdsEricson1.nex") 
# phylo <- ape::read.nexus("C:/git/STOC/Variables/data/THV/phylogeny/Ericson/AllBirdsEricson1.nex")
# phylo <- ape::read.nexus("C:/Users/Nathalie Adenot/STOC/cluster/data/AllBirdsEricson1.nex")

cat(" ------------------    Build consensus tree    -----------------------")

list_trees <- vector("list",100)
# select species from data for each tree  
for(i in 1:100){
  tree <- phylo[[i]]
  phylo_sp <- ape::drop.tip(tree, setdiff(tree$tip.label,sp_capt$phylo))
  list_trees[[i]] <- phylo_sp
}
ape::write.nexus(list_trees, file = "/scratchbeta/adenotn/data/test_subset_100.nex", translate = TRUE)
# ape::write.nexus(list_trees, file = "C:/git/STOC/Variables/data/THV/phylogeny/Ericson/test_subset_100.nex", translate = TRUE)
# ape::write.nexus(list_trees, file = "C:/Users/Nathalie Adenot/STOC/cluster/output/phylogeny/test_subset_1000.nex")

phylo_sub <- ape::read.nexus("/scratchbeta/adenotn/data/test_subset_100.nex")
# phylo_sub <- ape::read.nexus("C:/git/STOC/Variables/data/THV/phylogeny/Ericson/test_subset_100.nex")
# phylo_sub <- ape::read.nexus("C:/Users/Nathalie Adenot/STOC/cluster/output/phylogeny/test_subset_1000.nex")


# cons <- phangorn::consensusNet(phylo_sub)
# cons <- ape::consensus(phylo_sub, p = 1, check.labels = TRUE, rooted = TRUE)

# method 1 http://blog.phytools.org/2016/03/method-to-compute-consensus-edge.html 
t1<- phytools::consensus.edges(phylo_sub)
# plotTree(t1,fsize=0.4)


cat(" -----------    Phylogenetic variance-covariance matrix    ------------")
# A <- ape::vcv.phylo(cons) 
A <- ape::vcv.phylo(t1) 


# nom_especes <- fread("C:/git/STOC/Variables/data/THV/nom_especes.csv", encoding = "Latin-1")
nom_especes <- fread("/scratchbeta/adenotn/data/nom_especes.csv", encoding = "Latin-1")
# nom_especes <- fread("C:/Users/Nathalie Adenot/STOC/cluster/data/nom_especes.csv", encoding = "Latin-1")



cat(" -----------    Import THV    ------------")

# Import life-history traits
THV <- fread("/scratchbeta/adenotn/data/THV_all_sp.csv", encoding = "Latin-1")
# THV <- fread("C:/Users/Nathalie Adenot/STOC/cluster/data/THV_all_sp.csv", encoding = "Latin-1")
# THV <- fread("C:/git/STOC/cluster/data/THV_all_sp.csv", encoding = "Latin-1")

# transform species name
coef_sp_brms <- coef_sp %>% 
  left_join(nom_especes, by = c("Species" = "SP")) %>% 
  dplyr::select(Estimate, `Std. Error`, var, var_name, Species, nom_sc) %>% 
  mutate(nom_sc = gsub(" ", "_", nom_sc)) %>% 
  left_join(THV, by = c("Species" = "ESPECE")) %>% #add THV
  rename(Broods = `Broods per year`) %>% 
  mutate(nom_sc = gsub("Cyanistes", "Parus", nom_sc)) # in the tree, Cyanistes caeruleus is Parus caeruleus

colnames(coef_sp_brms)[1:7] <- c("est", "se", "var", "var_name", "ESPECE", "phylo", "nom_sc")

# check that all species are in the phylogeny  
diff <- setdiff(unique(coef_sp_brms$phylo), colnames(A))

# delete sp that are not in the tree from the df
coef_sp_brms <- coef_sp_brms %>% filter(!(phylo %in% diff))


# relevel factor  
coef_sp_brms <- coef_sp_brms %>% 
  mutate(HABITAT_SP = as.factor(HABITAT_SP)) %>% 
  mutate(MIGRATION = as.factor(MIGRATION))

coef_sp_brms$HABITAT_SP <- relevel(coef_sp_brms$HABITAT_SP, ref = "Terrestre")



cat(" ------------------    Loop on variables    -----------------------")
for(x in unique(coef_sp_brms$var_name)){
  cat(paste0("brms on ", x, "\n"))
  df <- subset(coef_sp_brms, var_name == x)
  df$obs <- 1:nrow(df)
  
  cat(paste0("MODEL WITH ALL TRAITS \n"))
  
  model_meta_all <- brm(est | se(se) ~ HABITAT_SP + MIGRATION + Broods + sti_europe +
                          (1|gr(phylo, cov = A)) + (1|obs),
                        data = df, family = gaussian(),
                        data2 = list(A = A),
                        prior = c(prior(normal(0, 10), "Intercept"),
                                  prior(student_t(3, 0, 10), "sd")),
                        control = list(adapt_delta = 0.95),
                        chains = 2, cores = 2, iter = 4000, warmup = 1000)
  
  saveRDS(model_meta_all, paste0("/scratchbeta/adenotn/output/model_meta_all_", x, ".rds"))
  
  print(summary(model_meta_all))
  
  png(paste0("/scratchbeta/adenotn/output/plot_model_meta_all_", x, ".png"))
  plot(model_meta_all)
  dev.off()
  
  
  
  cat(paste0("MODEL HABITAT_SP \n"))
  
  model_meta_hab <- brm(est | se(se) ~ HABITAT_SP +
                          (1|gr(phylo, cov = A)) + (1|obs),
                        data = df, family = gaussian(),
                        data2 = list(A = A),
                        prior = c(prior(normal(0, 10), "Intercept"),
                                  prior(student_t(3, 0, 10), "sd")),
                        control = list(adapt_delta = 0.95),
                        chains = 2, cores = 2, iter = 4000, warmup = 1000)
  
  saveRDS(model_meta_hab, paste0("/scratchbeta/adenotn/output/model_meta_hab", x, ".rds"))
  
  print(summary(model_meta_hab))
  
  png(paste0("/scratchbeta/adenotn/output/plot_model_meta_hab", x, ".png"))
  plot(model_meta_hab)
  dev.off()
  
  
  
  cat(paste0("MODEL MIGRATION \n"))
  
  model_meta_migr <- brm(est | se(se) ~ MIGRATION +
                           (1|gr(phylo, cov = A)) + (1|obs),
                         data = df, family = gaussian(),
                         data2 = list(A = A),
                         prior = c(prior(normal(0, 10), "Intercept"),
                                   prior(student_t(3, 0, 10), "sd")),
                         control = list(adapt_delta = 0.95),
                         chains = 2, cores = 2, iter = 4000, warmup = 1000)
  
  saveRDS(model_meta_migr, paste0("/scratchbeta/adenotn/output/model_meta_migr", x, ".rds"))
  
  print(summary(model_meta_migr))
  
  png(paste0("/scratchbeta/adenotn/output/plot_model_meta_migr", x, ".png"))
  plot(model_meta_migr)
  dev.off()
  
  cat(paste0("MODEL NB BROODS \n"))
  
  model_meta_broods <- brm(est | se(se) ~ Broods +
                             (1|gr(phylo, cov = A)) + (1|obs),
                           data = df, family = gaussian(),
                           data2 = list(A = A),
                           prior = c(prior(normal(0, 10), "Intercept"),
                                     prior(student_t(3, 0, 10), "sd")),
                           control = list(adapt_delta = 0.95),
                           chains = 2, cores = 2, iter = 4000, warmup = 1000)
  
  saveRDS(model_meta_broods, paste0("/scratchbeta/adenotn/output/model_meta_broods", x, ".rds"))
  
  print(summary(model_meta_broods))
  
  png(paste0("/scratchbeta/adenotn/output/plot_model_meta_broods", x, ".png"))
  plot(model_meta_broods)
  dev.off()
  
  
  
  cat(paste0("MODEL WITH STI \n"))
  
  model_meta_sti <- brm(est | se(se) ~ sti_europe +
                          (1|gr(phylo, cov = A)) + (1|obs),
                        data = df, family = gaussian(),
                        data2 = list(A = A),
                        prior = c(prior(normal(0, 10), "Intercept"),
                                  prior(student_t(3, 0, 10), "sd")),
                        control = list(adapt_delta = 0.95),
                        chains = 2, cores = 2, iter = 4000, warmup = 1000)
  
  saveRDS(model_meta_sti, paste0("/scratchbeta/adenotn/output/model_meta_sti_", x, ".rds"))
  
  print(summary(model_meta_sti))
  
  png(paste0("/scratchbeta/adenotn/output/plot_model_meta_sti_", x, ".png"))
  plot(model_meta_sti)
  dev.off()
  
  
  cat("DONE \n \n")
}