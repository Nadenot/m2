library(tidyverse)
library(ggplot2)
library(data.table)
library(lubridate)


## Data exploration
data_mod <- fread("C:/git/STOC/Variables/data/model/data_model.csv")



##########################    Spatial distribution    ########################## 
require(rgdal)
require(ggmap)
require(maptools)
library(sf)
require(maps)
# fond de carte
france <- map_data("france")
world1 <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))
france <- sf::st_as_sf(map('france', plot = FALSE, fill = TRUE))
mytheme <- theme(plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
                 panel.background = element_rect(colour = NA),
                 plot.background = element_rect(colour = NA),
                 axis.title = element_text(face = "bold",size = rel(1)),
                 axis.title.y = element_text(angle=90,vjust =2),
                 axis.title.x = element_text(vjust = -0.2),
                 legend.position=NULL)

## mean NDVI 
# All years
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord,aes(LON,LAT, colour=mNDVI),size=2)
gg <- gg + scale_color_gradient2(midpoint=0.6,  low="#f7fcb9", mid="#78c679", high="#004529", space = "Lab")
plot(gg)

# 2022
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord[data_mod_coord$YEAR == 2022,],aes(LON,LAT, colour=mNDVI),size=2)
plot(gg)

## mean SPEI 
# All years
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord,aes(LON,LAT, colour=mSPEI),size=2)
gg <- gg + scale_color_gradient2(midpoint=0,  low="red", mid="grey", high="blue", space = "Lab")
plot(gg)
# 2022
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord[data_mod_coord$YEAR == 2022,],aes(LON,LAT, colour=mSPEI),size=2)
gg <- gg + scale_color_gradient2(midpoint=0,  low="red", mid="grey", high="blue", space = "Lab")
plot(gg)

## sum SPEI 
# All years
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord,aes(LON,LAT, colour=sumSPEI),size=2)
gg <- gg + scale_color_gradient2(midpoint=0,  low="red", mid="grey", high="blue", space = "Lab")
plot(gg)
# 2022
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = data_mod_coord[data_mod_coord$YEAR == 2022,],aes(LON,LAT, colour=sumSPEI),size=2)
gg <- gg + scale_color_gradient2(midpoint=0,  low="red", mid="grey", high="blue", space = "Lab")
plot(gg)

## Productivity  
# Mean productivity per site (for all species)  
# All years
mean_prod <- data_mod %>% 
  group_by(ID_PROG) %>% 
  summarise(mean_prod = mean(Prod, na.rm = TRUE)) %>% 
  left_join(data_mod[,c("ID_PROG", "LON", "LAT")])
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = mean_prod,aes(LON,LAT, colour=mean_prod),size=2.5)
gg <- gg + scale_color_gradient2(midpoint=0.3,  low="#92c5de", mid="#ffe987", high="#d7191c", space = "Lab")
plot(gg)
# 2022
mean_prod <- data_mod %>% 
  filter(YEAR == 2022) %>% 
  group_by(ID_PROG) %>% 
  summarise(mean_prod = mean(Prod, na.rm = TRUE)) %>% 
  left_join(data_mod[,c("ID_PROG", "LON", "LAT")])
gg <- ggplot()
gg <- gg + geom_sf(data = world1,fill="white", colour="#7f7f7f", size=0.2)+ geom_sf(data = france,fill="white", colour="#7f7f7f", size=0.5)
gg <- gg + coord_sf(xlim=c(-5,9),ylim=c(41.5,52))
gg <- gg + geom_point(data = mean_prod,aes(LON,LAT, colour=mean_prod),size=2.5)
gg <- gg + scale_color_gradient2(midpoint=0.3,  low="#92c5de", mid="#ffe987", high="#d7191c", space = "Lab")
plot(gg)


# We want to know if there is a latitude effect on latitude  
# We don't really see something on the map, but I will test it:
m <- glmer(cbind(JUV,AD) ~ LAT + (1|ID_PROG) + (1|YEAR)+ (1|ESPECE), family = binomial, data = data_mod)
summary(m)


#####################     Temporal autocorrelation     #########################

## Plot evolution of NDVI  
gg <- ggplot(data_mod, aes(x = as.factor(YEAR), y = mNDVI)) +
  geom_boxplot() +
  scale_x_discrete(limits = factor(2000:2022), name = "Year") +
  theme(axis.text.x = element_text(angle=90)) 
plot(gg)

## Plot evolution of SPEI  
gg <- ggplot(data_mod, aes(x = as.factor(YEAR), y = sumSPEI)) +
  geom_boxplot() +
  scale_x_discrete(limits = factor(2000:2022), name = "Year") +
  theme(axis.text.x = element_text(angle=90)) 
plot(gg)

## Test for temporal autocorrelation


### Data exploration

library(ggplot2) # graph package
library(tinytex) # Pour la sortie pdf
library(corrplot)# Correlation matrix calculus
library(plot3D)# For 3D plot
library(DHARMa)# Model diagnosis
library(rcompanion)# Model pseudo R²
library(lattice)# multipanel graphics


################     check outliers and distribution      ######################
# Mean NDVI
par(mfrow=c(1,3))
# Length
# Cleveland plot
dotchart(data_mod$mNDVI,pch=16,col='blue',xlab='Mean spring NDVI') #OK
# Histogram
hist(data_mod$mNDVI,col='blue',xlab="Mean spring NDVI",main="") #not symmetrical
# Quantile-Quantile plot
qqnorm(data_mod$mNDVI,pch=16,col='blue',xlab='')
qqline(data_mod$mNDVI,col='red') #no normality

# SPEI
par(mfrow=c(1,3))
# Length
# Cleveland plot
dotchart(data_mod$sumSPEI,pch=16,col='blue',xlab='Total spring SPEI') #OK
# Histogram
hist(data_mod$sumSPEI,col='blue',xlab="Total spring SPEI",main="") #normal
# Quantile-Quantile plot
qqnorm(data_mod$sumSPEI,pch=16,col='blue',xlab='')
qqline(data_mod$sumSPEI,col='red') #normal


#############   Analysis of the potential relationships Y vs Xs   ##############

# For the graph: create a column prod = JUV/AD
data_mod_prod <- data_mod %>% 
  mutate(prod = JUV/(JUV+AD))
data_mod_prod$prod[is.infinite(data_mod_prod$prod)] <- NA


# NDVI
plot(data_mod_prod$prod~data_mod_prod$mNDVI,pch=16,col='blue',xlab='mean spring NDVI',ylab='Productivity (%JUV)')
# Add regression  
require(stats)
reg<-lm(prod ~ mNDVI, data = data_mod_prod)
coeff=coefficients(reg)
# Equation de la droite de regression :  
eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))
# Graphe
sp <- ggplot(data=data_mod_prod, aes(x=mNDVI, y=prod)) + geom_point(size = 0.7)
sp <- sp + geom_smooth(method="lm", se=FALSE)
sp <- sp + geom_abline(intercept = coeff[2], slope = coeff[1], color = "red", linewidth =1.5)+
  ggtitle(eq)
plot(sp)

# SPEI
plot(data_mod_prod$prod~data_mod_prod$sumSPEI,pch=16,col='blue',xlab='Total spring SPEI',ylab='Productivity (JUV/AD)')
# Add regression  
require(stats)
reg<-lm(prod ~ sumSPEI, data = data_mod_prod)
coeff=coefficients(reg)
# Equation de la droite de regression :  
eq = paste0("y = ", round(coeff[2],1), "*x + ", round(coeff[1],1))
# Graphe
sp <- ggplot(data=data_mod_prod, aes(x=sumSPEI, y=prod)) + geom_point(size = 0.7)
sp <- sp + geom_smooth(method="lm", se=FALSE)
sp <- sp + geom_abline(intercept = coeff[2], slope = coeff[1], color = "red", linewidth =1.5)+
  ggtitle(eq)
plot(sp)

###############         Check for colinearity between Xs       #################
data_mod <- fread("C:/git/STOC/Variables/data/model/data_model_all_a.csv")

library(DataExplorer)
plot_missing(data_mod[,9:23]) #plot % missing values for each variable
plot_correlation(data_mod[,9:23], cor_args = list("use" = "pairwise.complete.obs"), type="c") 
# correlation coefficient of 0.33 so OK


## Collinearity NDVI & SPEI
gg <- ggplot(data_mod, aes(y = mNDVI, x = sumSPEI)) +
  geom_point(stat = "identity")
plot(gg)              

gg <- ggplot(data_mod, aes(y = mNDVI, x = mSPEI)) +
  geom_point(stat = "identity")
plot(gg)              
# it looks thet same whether we use the sum or the mean

plot_missing(data_mod_prod) #plot % missing values for each variable
plot_correlation(data_mod_prod, cor_args = list("use" = "pairwise.complete.obs"), type="c") 


## Collinearity between mNDVI & aNDVI  
gg <- ggplot(data_mod, aes(y = yaNDVI, x = myNDVI)) +
  geom_point(stat = "identity")
plot(gg)              


## Collinearity between mTemp & aTemp  
gg <- ggplot(data_mod, aes(y = yaTemp, x = myTemp)) +
  geom_point(stat = "identity")
plot(gg)              


