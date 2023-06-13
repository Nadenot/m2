require(tidyverse)
require(data.table)
require(ggplot2)
require(FactoMineR)
require(Factoshiny)
require(factoextra)
require(RColorBrewer)

var_model_scaled <- fread("C:/git/STOC/cluster/data/var_model_final_allsp_scaled.csv")


### All variables  
## extraction var actives  
var_model_scaled_act <- var_model_scaled[, c(12:15, 17:18, 21:23, 26:27, 30:33)]
var_model_scaled_act <- var_model_scaled_act %>% 
  filter(!is.na(mean_aNDVI_early)) %>%
  filter(!is.na(early_SPEI)) %>% 
  filter(!is.na(early_nb_ECE_SPEI))

# PCA

res.pca <- prcomp(var_model_scaled_act, scale = FALSE)

summary(res.pca)
# screeplot(res.pca, las = 1, main = "", cex.lab = 1.5, xlab = "Principal components")
# 
# # biplot 
# biplot(res.pca, cex = 0.7)
# print(res.pca)
# 

# #eigenvalues 
# fviz_eig(res.pca)
# 
# 
# # graph individuals  
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Colorer par le cos2
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE)
# 
# # graph variables  
# fviz_pca_var(res.pca,
#              col.var = factor(c("Temperature", "Temperature", "SPEI", "SPEI", "NDVI", "NDVI", "Temperature", "NDVI", "SPEI", "Temperature", "SPEI", "Temperature", "Temperature", "SPEI", "SPEI")),
#              # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     
# )
# 
# 
# # contribution des variables aux axes
# var <- get_pca_var(res.pca)
# head(var$contrib, 4)
# # Contributions des variables à PC1
# fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
# # Contributions des variables à PC2
# fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
# 


## Colorer par groupe 
var_model_scaled_pca <- var_model_scaled[, c(1:2, 12:15, 17:18, 21:23, 26:27, 30:33)] %>% distinct() %>% 
  filter(!is.na(mean_aNDVI_early)) %>%
  filter(!is.na(early_SPEI)) %>% 
  filter(!is.na(early_nb_ECE_SPEI)) %>% 
  mutate(YEAR = as.factor(YEAR)) %>% 
  rename(aNDVI_early = mean_aNDVI_early, aNDVI_late = mean_aNDVI_late, aT_early = mean_aT_early, aT_late = mean_aT_late, 
         SPEI_early = early_SPEI, SPEI_late =late_SPEI, ECE_SPEI = nb_ECE_SPEI, ECE_aT = nb_ECE_aTemp, 
         ECE_SPEI_early = early_nb_ECE_SPEI, ECE_SPEI_late = late_nb_ECE_SPEI, ECE_aT_early = early_nb_ECE_aTemp, ECE_aT_late = late_nb_ECE_aTemp)

var_model_scaled_pca_y <- var_model_scaled_pca %>% filter(YEAR %in% c("2003", "2020", "2021", "2022"))

iris.pca <- PCA(var_model_scaled_pca_y[, - c(1:2)], graph = FALSE)

# fviz_pca_ind(iris.pca,
#              geom.ind = "point", # Montre les points seulement (mais pas le "text")
#              col.ind = var_model_scaled_pca_y$YEAR, # colorer by groups
#              palette = "lancet",
#              addEllipses = TRUE, # Ellipses de concentration
#              legend.title = "Groups"
# )
# 
# # Ajoutez des ellipses de confiance
# fviz_pca_ind(iris.pca, geom.ind = "point", col.ind = var_model_scaled_pca_y$YEAR, 
#              # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
#              addEllipses = TRUE, ellipse.type = "confidence",
#              legend.title = "Groups"
# )
# 
# # barycentres
# fviz_pca_ind (iris.pca,
#               geom.ind = "point", # afficher les points seulement (pas de "texte")
#               col.ind = var_model_scaled_pca_y$YEAR, # Couleur par groupes
#               legend.title = "Groupes",
#               mean.point = FALSE)


# biplot by group
fviz_pca_biplot (iris.pca,
                 col.ind = var_model_scaled_pca_y$YEAR, #palette = "simpsons",
                 addEllipses = TRUE, label = "var",
                 col.var = factor(c("Temperature", "Temperature", "SPEI", "SPEI", "NDVI", "NDVI", "Temperature", "NDVI", "SPEI", "Temperature", "SPEI", "Temperature", "Temperature", "SPEI", "SPEI")),
                 repel = TRUE,
                 legend.title = list(fill = "Year", color = "Clusters"))+
ggpubr::fill_palette(c("#8dd3c7", "#fdc086", "#bebada", "#fbb4ae"))+      # Couleur des individus
ggpubr::color_palette(c("#8dd3c7", "#fdc086", "#bebada", "#fbb4ae","#7fc97f", "#1f78b4", "#ec7063"))      # Couleur des variables

# 
# fviz_pca_biplot(iris.pca,
#                 # Colueur de remplissage des individdus par groupes
#                 geom.ind = "point",
#                 pointshape = 21,
#                 pointsize = 2.5,
#                 fill.ind = var_model_scaled_pca_y$YEAR,
#                 col.ind = "black",
#                 # Colorer les variables par groupes
#                 col.var = factor(c("Temperature", "Temperature", "SPEI", "SPEI", "NDVI", "NDVI", "Temperature", "NDVI", "SPEI", "Temperature", "SPEI", "Temperature", "Temperature", "SPEI", "SPEI")),
#                 legend.title = list(fill = "Year", color = "Clusters"),
#                 repel = TRUE        # Evite le chévauchement du texte
# )#+
#   # ggpubr::fill_palette("jco")+      # Couleur des individus
#   # ggpubr::color_palette("aaas")      # Couleur des variables


## Supplementary variables  
# groups <- as.factor(var_model_scaled_pca_y[,c(1:3, 7:8)])
# fviz_pca_ind(res.pca,
#              col.ind = groups, # colorer par groupes
#              palette = c("#00AFBB",  "#FC4E07"),
#              addEllipses = TRUE, # Ellipse de concentration
#              ellipse.type = "confidence",
#              legend.title = "Groups",
#              repel = TRUE
# )

## Individus et variables supplémentaires 
var_model_scaled_supp <- var_model_scaled[, c(2, 9:11, 12:15, 16:18, 21:23, 26:27, 30:33)] %>% distinct() %>% 
  filter(!is.na(mean_aNDVI_early)) %>%
  filter(!is.na(early_SPEI)) %>% 
  filter(!is.na(early_nb_ECE_SPEI)) %>% 
  mutate(YEAR = as.factor(YEAR)) %>% 
  rename(aNDVI_early = mean_aNDVI_early, aNDVI_late = mean_aNDVI_late, aT_early = mean_aT_early, aT_late = mean_aT_late, 
         SPEI_early = early_SPEI, SPEI_late =late_SPEI, ECE_SPEI = nb_ECE_SPEI, ECE_aT = nb_ECE_aTemp, 
         ECE_SPEI_early = early_nb_ECE_SPEI, ECE_SPEI_late = late_nb_ECE_SPEI, ECE_aT_early = early_nb_ECE_aTemp, ECE_aT_late = late_nb_ECE_aTemp)

  
res.pca <- PCA(var_model_scaled_supp,  
               quanti.sup = c(2:3, 9), quali.sup = c(1,4), graph=FALSE)
# quanti
# fviz_pca_var(res.pca)
# Changer la couleur des variables
fviz_pca_var(res.pca,
             col.var = factor(c("Temperature", "Temperature", "SPEI", "SPEI", "NDVI", "NDVI", "Temperature", "NDVI", "SPEI", "Temperature", "SPEI", "Temperature", "Temperature", "SPEI", "SPEI")), # Variables actives
             col.quanti = "black", # variables quantitatives supl.
             repel = TRUE
)+
  ggpubr::color_palette(c("#7fc97f", "#1f78b4", "#ec7063"))      # Couleur des variables
#ggpubr::color_palette(c("#ec7063","#ec7063", "#1f78b4", "#1f78b4","#7fc97f","#7fc97f","#ec7063", "#7fc97f", "#1f78b4", "#ec7063", "#1f78b4","#ec7063", "#ec7063", "#1f78b4","#1f78b4"))      # Couleur des variables

# # Cacher les variables actives sur le graphique,
# # ne montrent que des variables supplémentaires
# fviz_pca_var(res.pca, invisible = "var")
# # Cacher les variables supplémentaires
# fviz_pca_var(res.pca, invisible = "quanti.sup")

###################### Site characteristics  ###############################
# 
# ## extraction var actives  
# var_model_scaled_act <- var_model_scaled[, c(9:10, 16)]
# var_model_scaled_act <- var_model_scaled_act %>% 
#   filter(!is.na(meanNDVI)) 
# 
# # Factoshiny(var_model_scaled_act)
# # PCA
# 
# res.pca <- prcomp(var_model_scaled_act, scale = FALSE)
# 
# summary(res.pca)
# # screeplot(res.pca, las = 1, main = "", cex.lab = 1.5, xlab = "Principal components")
# # 
# # # biplot 
# # biplot(res.pca, cex = 0.7)
# # print(res.pca)
# 
# 
# #eigenvalues 
# fviz_eig(res.pca)
# 
# 
# # graph individuals  
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Colorer par le cos2
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE)
# 
# # graph variables  
# fviz_pca_var(res.pca,
#              col.var = "contrib", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     
# )
# 
# 
# ## Supplementary variables  
# groups <- var_model_scaled %>%
#   filter(!is.na(meanNDVI)) %>% 
#   select(c(1,3, 7:8)) %>% 
#   mutate_all(as.factor)
# fviz_pca_ind(res.pca,
#              col.ind = groups, # colorer par groupes
#              palette = c("#00AFBB",  "#FC4E07", "#E7B800", "#252525"),
#              addEllipses = TRUE, # Ellipse de concentration
#              ellipse.type = "confidence",
#              legend.title = "Groups",
#              repel = TRUE
# )
# 
# 
# 
# 
# 
# ######################  Only spring  #########################
# 
# ## extraction var actives  
# var_model_scaled_act <- var_model_scaled[, c(9:10, 12:18, 21:23, 26:33)]
# var_model_scaled_act <- var_model_scaled_act %>% 
#   filter(!is.na(meanNDVI)) %>%
#   filter(!is.na(early_SPEI)) %>% 
#   filter(!is.na(early_nb_ECE_SPEI))
# 
# # PCA
# 
# res.pca <- prcomp(var_model_scaled_act, scale = FALSE)
# 
# summary(res.pca)
# screeplot(res.pca, las = 1, main = "", cex.lab = 1.5, xlab = "Principal components")
# 
# # biplot 
# biplot(res.pca, cex = 0.7)
# print(res.pca)
# 
# 
# #eigenvalues 
# fviz_eig(res.pca)
# 
# 
# # graph individuals  
# fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Colorer par le cos2
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE)
# 
# # graph variables  
# fviz_pca_var(res.pca,
#              col.var = "contrib", 
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     
# )
# 
# 
# ## Supplementary variables  
# groups <- as.factor(var_model_scaled[,c(1:3, 7:8)])
# fviz_pca_ind(res.pca,
#              col.ind = groups, # colorer par groupes
#              palette = c("#00AFBB",  "#FC4E07"),
#              addEllipses = TRUE, # Ellipse de concentration
#              ellipse.type = "confidence",
#              legend.title = "Groups",
#              repel = TRUE
# )