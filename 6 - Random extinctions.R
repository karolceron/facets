#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires

## Loading R packages

library("BAT") 
library(ape)
library(vegan)
library(ade4)
library(geiger)
library(picante)
library(phytools)
library(geometry)
library(hypervolume)
library("bipartite")
library("iNEXT")
library("vegan")
library("ggplot2")
library('igraph')
library("mFD")
library(rgl)
library(dplyr)

############################
## null model - random extinctions
###########################

temp <- dir("./data/nets")
paths <- paste0("./data/nets/",temp)
myfiles = lapply(paths, read.table)
names(myfiles) <- temp

int <- t(myfiles[[4]]) #choose ecoregion 1 to 4
int <- as.data.frame(int)
rowSums(int) #checking data
colSums(int) #checking data


int <- int[, which(colSums(int) != 0)] 
int 
 
#transforming data 
int$Formicidae <- log(1+int$Formicidae) 
int$Isoptera <- log(1+int$Isoptera) 
int

#creating 100 lists 
ext.ran <- list()
for (i in 1:100) {
  ext.ran[[i]]<-int}

data <- data.frame(matrix(ncol = 3, nrow = 20))
a <- c("int_vol", "fun_vol", "phy_vol")
colnames(data) <- a
data

vol_int <- list()
vol_fun <- list()
vol_phy <- list()

ext <- 0

### functions
hullvol <- function(x) {  x$vol
}

teste <- function(x) {
  if (ncol(x$vectors)>2){
    x$vectors[,1:3] 
  }}

teste2 <- function(x) {
  if (ncol(x$vectors)>2){
    x$vectors[,1:2] 
  }}

trai <- function(x) {
  traits[which(rownames(traits) %in% rownames(x)), ]  
}

zero <- function(x) {
  x[, which(colSums(x) != 0)]  
}

um <- function(x) {
  x[, which(colSums(x) != nrow(x))]  
}

dis <- function(x) {
  dist.ktab(x, type= c("Q", "Q", "Q","B", "O"))}

###########

for (i in 1:nrow(data)) {
  
  ext <- ext+1
  ext.ran <- lapply(ext.ran, extinction, participant="lower", method="random")
  ext.ran <- lapply(ext.ran, empty)
  
  
  #interaction space
  int <- lapply(ext.ran, vegdist, method="horn")
  pcoa <- lapply(int, pcoa)
  
  
  #### interaction space
  vector <- lapply(pcoa,teste)
  vector <- vector[lengths(vector)!=0]
  hull <- lapply (vector,convhulln,options="Fa", output.options=TRUE) 
  vol <- lapply (hull,hullvol) 
  mean_int <- mean(unlist(vol))
  data$int_vol[ext] <- mean_int
  data$int_vol
  
  df <- as.data.frame(unlist(vol))
  df
  vol_int[ext] <- df
  vol_int
  
  
  #functional space
  traits <- read.table("./data/traits.txt", header = T)
  traits <- as.data.frame(traits)
  trait_r <- list()
  trait_r <- lapply(ext.ran,trai) 
  trait_r <- lapply(trait_r, zero)
  
  trait_r_habitat=list()
  for (i in 1:100){trait_r_habitat[[i]]=trait_r[[i]][,-which(colnames(trait_r[[i]])%in% c("size","clutch_size","mass","mode"))]
  }
  
  
  habitat <- function(x) {
    prep.binary(x,col.blocks=ncol(x),label="habitat")}
  
  
  hab <- lapply(trait_r_habitat,habitat)
  
  massa <- function(x) {
    data.frame(body=x$mass)}
  mass<-lapply(trait_r, massa)
  
  clu <- function(x) {
    data.frame(clutch=x$clutch_size)}
  clutch<-lapply(trait_r, clu)
  
  siz <- function(x) {
    data.frame(size=x$size)}
  size<-lapply(trait_r, siz)
  
  mod <- function(x) {
    data.frame(mode=x$mode)}
  mode<-lapply(trait_r, mod)
  
  ktab <- list()
  for (i in 1:100) {
    ktab[[i]] <-ktab.list.df(list(size[[i]], clutch[[i]], mass[[i]],hab[[i]],mode[[i]]), rownames = rownames(trait_r[[i]]))}
  
  
  int <- lapply(ktab, dis)
  
  int2 <- int[sapply(int, function(x) all(!is.na(x)))]
  
  
  pcoa <- lapply(int2, pcoa)
  vector <- lapply(pcoa,teste)
  vector <- vector[lengths(vector)!=0]
  hull <- lapply (vector,convhulln,options="Fa", output.options=TRUE) 
  vol <- lapply (hull,hullvol) 
  mean_fun <- mean(unlist(vol))
  data$fun_vol[ext] <- mean_fun
  data$fun_vol
  
  df <- as.data.frame(unlist(vol))
  df
  vol_fun[ext] <- df
  
  
  #phylogenetic space
  tree <- read.tree("./data/tree.nwk")
  tree
  
  name.check(tree, trait_r[[1]])
  
  prunedphy <- prune.sample(t(trait_r[[1]]),tree)
  
  for (i in 1:length(trait_r)) {
    prunedphy[[i]] <- prune.sample(t(trait_r[[i]]),tree)}
  
  cophe <- function(x) {
    cophenetic(x)}
  
  dist <- lapply(prunedphy, cophe)
  pcoa <- lapply(dist, pcoa)
  vector <- lapply(pcoa,teste)
  vector <- vector[lengths(vector)!=0] 
  hull <- lapply (vector,convhulln,options="Fa", output.options=TRUE) 
  vol <- lapply (hull,hullvol) 
  mean_phy <- mean(unlist(vol))
  data$phy_vol[ext] <- mean_phy
  data$phy_vol
  
  df <- as.data.frame(unlist(vol))
  df
  vol_phy[ext]<-df
  
}


write.csv(data, "random.csv")

#########################################
##### plotting diversities in random extinctions
######################################

data <- read.table("./data/desv.txt", header=T) #diversities after random extinctions for all ecoregions
data


### run the code below to each diversity facet
#### example with interaction diversity

a <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(type %in% "random") %>%
  filter(area %in% "cerrado") %>%
  ggplot() +
  aes(x = ext, y = mean, na.rm = TRUE) +
  geom_ribbon(aes(ymin=desv1,ymax=desv2),alpha=0.3)+
  geom_line(size = 0.5, linetype ="dashed", color="white")+
  geom_point(size = 2,color="white")+
  theme_classic(base_size = 18) +
  ylim(0, 0.3)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Interaction diversity",
  )  +scale_colour_manual("",values="burlywoord3")
a

land <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(area %in% "cerrado") %>%
  filter(type %in% c("land"))

cerrado <- a + geom_line(y = land$mean, size = 1, na.rm = TRUE)


b <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(type %in% "random") %>%
  filter(area %in% "chaco") %>%
  ggplot() +
  aes(x = ext, y = mean, na.rm = TRUE) +
  geom_ribbon(aes(ymin=desv1,ymax=desv2),alpha=0.3)+
  geom_line(size = 0.5, linetype ="dashed", color="white")+
  geom_point(size = 2,color="white")+
  theme_classic(base_size = 18) +
  ylim(0, 0.3)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Interaction diversity",
  )  +scale_colour_manual("",values="burlywoord3")
b

land <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(area %in% "chaco") %>%
  filter(type %in% c("land"))

chaco <- b + geom_line(y = land$mean, size = 1, na.rm = TRUE)


c <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(type %in% "random") %>%
  filter(area %in% "ma") %>%
  ggplot() +
  aes(x = ext, y = mean, na.rm = TRUE) +
  geom_ribbon(aes(ymin=desv1,ymax=desv2),alpha=0.3)+
  geom_line(size = 0.5, linetype ="dashed", color="white")+
  geom_point(size = 2,color="white")+
  theme_classic(base_size = 18) +
  ylim(0, 0.3)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Interaction diversity",
  )  +scale_colour_manual("",values="burlywoord3")
c

land <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(area %in% "ma") %>%
  filter(type %in% c("land"))

ma <- c + geom_line(y = land$mean, size = 1, na.rm = TRUE)

d <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(type %in% "random") %>%
  filter(area %in% "pantanal") %>%
  ggplot() +
  aes(x = ext, y = mean, na.rm = TRUE) +
  geom_ribbon(aes(ymin=desv1,ymax=desv2),alpha=0.3)+
  geom_line(size = 0.5, linetype ="dashed", color="white")+
  geom_point(size = 2,color="white")+
  theme_classic(base_size = 18) +
  ylim(0, 0.3)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Interaction diversity",
  )  +scale_colour_manual("",values="burlywoord3")
d

land <- data %>%
  filter(diversity %in% "interaction") %>% #change facet
  filter(area %in% "pantanal") %>%
  filter(type %in% c("land"))

pantanal <- d + geom_line(y = land$mean, size = 1, na.rm = TRUE)


grid.arrange(ma,cerrado, chaco,pantanal) # diversities to random extinctions
