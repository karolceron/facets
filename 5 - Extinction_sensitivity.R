#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires



##############################
##### sensitivity analysis
#############################


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
library(tnet)
library(dplyr)

#################################################################
## Loading the databases and preparing the data


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
int


int_vol <- data.frame(matrix(ncol = 6, nrow = nrow(int)))
fun_vol <- data.frame(matrix(ncol = 6, nrow = nrow(int)))
phy_vol <- data.frame(matrix(ncol = 6, nrow = nrow(int)))


vol_int <- list()
vol_fun <- list()
vol_phy <- list()


##### first round
ext <- 1
ext.suarede <- int

## Converting interaction to a distance matrix (horn distance) and transforming them with a PCoA
horn <- vegdist(ext.suarede, method="horn")
pcoa <- pcoa(horn)

#saving vectors
vec.list <-  list()
for(k in 2:6){
  pcoa_matrix <- pcoa
  n_vectors <- ncol(pcoa_matrix$vectors)
  num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
  vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
}


for(k in 2:6){
  c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
  vol <- c$vol #volume
  int_vol[ext,k] <- vol 
}


######################################
# functional space anurans 
########################################

#creating community matrix
anura <- row.names(ext.suarede)
site <- rep(1, nrow(ext.suarede)) 
comm <- data.frame(matrix(nrow = 1,data = site))
colnames(comm) <- anura
rownames(comm) <- 'site'
comm

traits <- read.table("data/traits.txt", header = T)
traits <- as.data.frame(traits)

traits <- traits[which(rownames(traits) %in% colnames(comm)), ] 
traits
dim(traits)

traits <- traits[, which(colSums(traits) != 0)] 
traits

## Converting traits to a distance matrix (gower distance) and transforming them with a PCoA
hab <- prep.binary(traits[,4:7],col.blocks=4,label="habitat")
mass <- data.frame(body=traits$mass)
clutch <- data.frame(clutch=traits$clutch_size)
size <- data.frame(size=traits$size)
mode <- data.frame(mode=traits$mode)
ktab <- ktab.list.df(list(size, clutch, mass,hab, mode), rownames = rownames(traits))
dist <- dist.ktab(ktab, type= c("Q", "Q", "Q","B", "O")) #gower mixed
pcoa <- pcoa(dist)



#saving vectors
vec.list <-  list()
for(k in 2:6){
  pcoa_matrix <- pcoa
  n_vectors <- ncol(pcoa_matrix$vectors)
  num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
  vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
}


for(k in 2:6){
  c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
  vol <- c$vol #volume
  fun_vol[ext,k] <- vol 
}


#################################
# phylogenetic space anurans 
###############################
tree <- read.tree("data/tree.nwk")
tree
comu <- t(comm)
dim(comu)
name.check(tree, comu) #checar nomes
prunedphy <-prune.sample(t(comu),tree) #cortar arvore
prunedphy
plotTree(prunedphy)


## Converting tree to a distance matrix (cophenetic distance) and transforming them with a PCoA
phydist <- cophenetic(prunedphy)
pcoa <- pcoa(phydist)  

vec.list <-  list()
for(k in 2:6){
  pcoa_matrix <- pcoa
  n_vectors <- ncol(pcoa_matrix$vectors)
  num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
  vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
}

for(k in 2:6){
  #phylogenetic space
  c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
  vol <- c$vol #volume
  phy_vol[ext,k] <- vol 
}


#==============================
#Extinction simulations
#==============================

# extinction order to 8.5 land use

ext_order <- read.table("./data/ext_order.txt", header = T)


####### choose extinction order according to the ecoregion #########


#ext_order <- ext_order %>%  filter(site %in% "cerrado") #ecoregion 1

#ext_order <- ext_order %>%  filter(site %in% "chaco") #ecoregion 2

#ext_order <- ext_order %>%  filter(site %in% "af") #ecoregion 3

ext_order <- ext_order %>%  filter(site %in% "pantanal") #ecoregion 4



for (i in 1:length(ext_order$species)) {
  
  ext<-ext+1
  ext.suarede <- ext.suarede[!(row.names(ext.suarede) %in% ext_order[ext,1]),] #removing in extinction order
  
  ## Converting interaction to a distance matrix (horn distance) and transforming them with a PCoA
  horn=vegdist(ext.suarede, method="horn")
  pcoa=pcoa(horn)
  
  
  vec.list <-  list()
  for(k in 2:6){
    pcoa_matrix <- pcoa
    n_vectors <- ncol(pcoa_matrix$vectors)
    num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
    vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
  }
  
  for(k in 2:6){
    c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
    vol <- c$vol #volume
    int_vol[ext,k] <- vol 
  }
  
  
  ######################################
  # functional space anurans traits
  ########################################
  
  #creating community matrix
  anura<-row.names(ext.suarede)
  site<-rep(1, nrow(ext.suarede)) 
  comm <- data.frame(matrix(nrow = 1,data = site))
  colnames(comm) <- anura
  rownames(comm) <- 'site'
  comm
  
  traits <- read.table("data/traits.txt", header = T)
  traits<-as.data.frame(traits)
  
  
  traits <- traits[which(rownames(traits) %in% colnames(comm)), ] 
  traits
  dim(traits)
  
  traits<-traits[, which(colSums(traits) != 0)] 
  traits
  

  ## Converting traits to a distance matrix (gower distance) and transforming them with a PCoA
  hab=prep.binary(traits[,4:7],col.blocks=4,label="habitat")
  mass=data.frame(body=traits$mass)
  clutch=data.frame(clutch=traits$clutch_size)
  size=data.frame(size=traits$size)
  mode=data.frame(mode=traits$mode)
  ktab=ktab.list.df(list(size, clutch, mass,hab, mode), rownames = rownames(traits))
  dist=dist.ktab(ktab, type= c("Q", "Q", "Q","B", "O")) #gower mixed
  pcoa=pcoa(dist)
 
  vec.list <-  list()
  for(k in 2:6){
    pcoa_matrix <- pcoa
    n_vectors <- ncol(pcoa_matrix$vectors)
    num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
    vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
  }
  
  for(k in 2:6){
    c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
    vol <- c$vol #volume
    fun_vol[ext,k] <- vol 
  }
  
  #################################
  # phylogenetic space anurans 
  ###############################
  tree=read.tree("data/tree.nwk")
  tree
  comu<-t(comm)
  dim(comu)
  name.check(tree, comu) #checar nomes
  prunedphy <-prune.sample(t(comu),tree) #cortar arvore
  prunedphy
  plotTree(prunedphy)
  
  
  ## Converting tree to a distance matrix (cophenetic distance) and transforming them with a PCoA
  phydist=cophenetic(prunedphy)
  pcoa=pcoa(phydist)  
  
  vec.list <-  list()
  for(k in 2:6){
    pcoa_matrix <- pcoa
    n_vectors <- ncol(pcoa_matrix$vectors)
    num_to_save <- ifelse(n_vectors >= k, k, n_vectors)
    vec.list[[k]] <- pcoa_matrix$vectors[, 1:num_to_save]
  }
  
  for(k in 2:6){
    c <- convhulln(vec.list[[k]],options="Fa", output.options=TRUE)
    vol <- c$vol #volume
    phy_vol[ext,k] <- vol 
  }
  
}

#formatting data

int_vol <- int_vol[,-1]
fun_vol <- fun_vol[,-1]
phy_vol <- phy_vol[,-1]

colnames(int_vol)<-c("int_vol2", "int_vol3","int_vol4","int_vol5","int_vol6")
colnames(fun_vol)<-c("fun_vol2", "fun_vol3","fun_vol4","fun_vol5","fun_vol6")
colnames(phy_vol)<-c("phy_vol2", "phy_vol3","phy_vol4","phy_vol5","phy_vol6")


data <- cbind(int_vol, fun_vol, phy_vol)
data


write.csv(data, "sensitivity.csv")



############################################
############# plotting sensitivity
#########################################

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(GGally)
library(viridis)
library(lattice)
library(cowplot)
library(wesanderson)
library(gridExtra)
library(vegan)

#loading sensitivity data from all sites
sensi <- read.table("./data/sensi.txt", header = T) 
sensi
head(sensi)


#standardizing values
sensi3 <- sensi %>%
  filter(axis %in% "three") 
sensi3[,3]<-decostand(sensi3[,3],method = 'range', na.rm = T)
sensi3

sensi4 <- sensi %>%
  filter(axis %in% "four") 
sensi4[,3]<-decostand(sensi4[,3],method = 'range', na.rm = T)
sensi4

sensi5 <- sensi %>%
  filter(axis %in% "five") %>%
  filter(ext %in% 1:12) 
sensi5[,3]<-decostand(sensi5[,3],method = 'range', na.rm = T)
sensi5

sensi6 <- sensi %>%
  filter(axis %in% "six") %>%
  filter(ext %in% 1:9) 
sensi6[,3]<-decostand(sensi6[,3],method = 'range', na.rm = T)
sensi6

#joining data
new <- rbind(sensi3, sensi4, sensi5, sensi6)
new

#plotting by retained axis
three <- new %>%
  #filter(site %in% "pantanal") %>%
  filter(axis %in% "three") %>%
  ggplot() +
  aes(x = ext, y = value, colour = type, linetype=site) +
  geom_line(size = 0.8)+
  geom_point(size = 1)+
  theme_classic(base_size = 18) +
  #ylim(0, 0.001)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Diversity",
  ) 

four <- new %>%
  #filter(site %in% "pantanal") %>%
  filter(axis %in% "four") %>%
  filter(ext %in% 1:12) %>%
  ggplot() +
  aes(x = ext, y = value, colour = type, linetype=site) +
  geom_line(size = 0.8)+
  geom_point(size = 1)+
  theme_classic(base_size = 18) +
  #ylim(0, 0.001)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Diversity",
  ) 


five <- new %>%
  #filter(site %in% "pantanal") %>%
  filter(axis %in% "five") %>%
  filter(ext %in% 1:12) %>%
  ggplot() +
  aes(x = ext, y = value, colour = type, linetype=site) +
  geom_line(size = 0.8)+
  geom_point(size = 1)+
  theme_classic(base_size = 18) +
  ylim(0, 1)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Diversity",
  ) 


six <- new %>%
  #filter(site %in% "pantanal") %>%
  filter(axis %in% "six") %>%
  filter(ext %in% 1:9) %>%
  ggplot() +
  aes(x = ext, y = value, colour = type, linetype=site) +
  geom_line(size = 0.8)+
  geom_point(size = 1)+
  theme_classic(base_size = 18) +
  #ylim(0, 0.001)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Diversity",
  ) 


plot_grid(three, four, five, six, labels=c("A)", "B)", "C)", "D"))


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

#creating 100 lists 
ext.ram <- list()
for (i in 1:100) {
  ext.ram[[i]]<-int}

data <- data.frame(matrix(ncol = 4, nrow = 20))
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
  ext.ram <- lapply(ext.ram, extinction, participant="lower", method="random")
  ext.ram <- lapply(ext.ram, empty)
  
  
  #interaction space
  int <- lapply(ext.ram, vegdist, method="horn")
  pcoa <- lapply(int, pcoa)
  
  
  #### interaction space
  vector <- lapply(pcoa,teste)
  vector <- vector[lengths(vector)!=0]
  hull <- lapply (vector,convhulln,options="Fa", output.options=TRUE) 
  vol <- lapply (hull,hullvol) 
  mean_int <- mean(unlist(vol))
  data$int_vol[ext] <- mean_fun
  data$int_vol
  
  df <- as.data.frame(unlist(vol))
  df
  vol_int[ext] <- df
  vol_int
  
  
  #functional space
  traits <- read.table("./data/traits.txt", header = T)
  traits <- as.data.frame(traits)
  trait_r <- list()
  trait_r <- lapply(ext.ram,trai) 
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

data <- read.table("./data/desv.txt", header=T) #diversities data to random extinctions
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
