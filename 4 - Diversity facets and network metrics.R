#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires

#loading packages
library("BAT") 
library("ape")
library("vegan")
library("ade4")
library("geiger")
library("picante")
library("phytools")
library("geometry")
library("bipartite")
library("tnet")
library("ggplot2")
library("wesanderson")
library("gridExtra")
library("dplyr")
library("cowplot")
library("mFD")

#read data

temp <- dir("./data/nets")
paths <- paste0("./data/nets/",temp)
myfiles = lapply(paths, read.table)
names(myfiles) <- temp

int <- t(myfiles[[4]]) #choose ecoregion 1 to 4
int <- as.data.frame(int)
rowSums(int) #checking data
colSums(int) #checking data

int <- int[, which(colSums(int) != 0)] #removing zero columns
int


int$Formicidae <- log(1+int$Formicidae) #transforming more abundant data
int$Isoptera <- log(1+int$Isoptera) #transforming more abundant data

######################################
# interaction diversity 
########################################

#creating community matrix
anura <- row.names(int)
site <- rep(1, nrow(int)) 
comm <- data.frame(matrix(nrow = 1,data = site))
colnames(comm) <- anura
rownames(comm) <- 'site'
comm

## Converting interaction to a distance matrix (horn distance) and transforming them with a PCoA
horn <- vegdist(int, method="horn")
pcoa <- pcoa(horn)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector <- pcoa$vectors[,1:3]

#convex hull 
c <- convhulln(vector,options="Fa", output.options=TRUE)
vol <- c$vol #volume
vol #interaction diversity/volume

# relation among traits and axis 1-10
tr_faxes3 <- mFD::traits.faxes.cor(
  sp_tr          = int[,1:10], 
  sp_faxes_coord = vector[ , c("Axis.1", "Axis.2")], 
  plot           = TRUE)

tr_faxes3$"tr_faxes_plot"

#only significant results
tr_faxes3$"tr_faxes_stat"[which(tr_faxes3$"tr_faxes_stat"$"p.value" < 0.05), ]


#plotting space
fspaces <- mFD::quality.fspaces(
  sp_dist             = horn,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

sp_faxes <- fspaces$"details_fspaces"$"sp_pc_coord"

big_plot1 <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes[ , c("PC1", "PC2", "PC3")],
  faxes           = c("PC1","PC2"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "black",
  fill_pool       = "white",
  shape_pool      = 3,
  size_pool       = 2,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "black",
  fill_vert       = "#7fbf7b",
  shape_vert      = 23,
  size_vert       = 3,
  plot_sp_nm      = rownames(vector),
  nm_size         = 2,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


big_plot1


######################################
# functional diversity 
########################################

traits <- read.table("data/traits.txt", header = T) #read data
traits <- as.data.frame(traits)

traits <- traits[which(rownames(traits) %in% colnames(comm)), ] #cutting traits according community 
traits
dim(traits)

traits <- traits[, which(colSums(traits) != 0)] #removing zero columns
traits

## Converting traits to a distance matrix (gower distance) and transforming them with a PCoA
hab <- prep.binary(traits[,4:6],col.blocks=3,label="habitat")
mass <- data.frame(body=traits$mass)
clutch <- data.frame(clutch=traits$clutch_size)
size <- data.frame(size=traits$size)
mode <- data.frame(mode=traits$mode)
ktab <- ktab.list.df(list(size, clutch, mass,hab, mode), rownames = rownames(traits))
dist <- dist.ktab(ktab, type= c("Q", "Q", "Q","B", "O")) #gower mixed
pcoa <- pcoa(dist)
barplot(pcoa$values$Eigenvalues)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector <- pcoa$vectors[,1:3]

#convex hull 
c <- convhulln(vector,options="Fa", output.options=TRUE)
vol <- c$vol
vol #functional diversity/volume


# relation among traits and axis
tr_faxes3 <- mFD::traits.faxes.cor(
  sp_tr          = traits, 
  sp_faxes_coord = vector[ , c("Axis.1", "Axis.2")], 
  plot           = TRUE)

tr_faxes3$"tr_faxes_plot"

#only significant results
tr_faxes3$"tr_faxes_stat"[which(tr_faxes3$"tr_faxes_stat"$"p.value" < 0.05), ]


fspaces <- mFD::quality.fspaces(
  sp_dist             = dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

sp_faxes <- fspaces$"details_fspaces"$"sp_pc_coord"

big_plot1 <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes[ , c("PC1", "PC2", "PC3")],
  faxes           = c("PC1","PC3"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "black",
  fill_pool       = "white",
  shape_pool      = 3,
  size_pool       = 2,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "black",
  fill_vert       = "#7fbf7b",
  shape_vert      = 23,
  size_vert       = 3,
  plot_sp_nm      = rownames(vector),
  nm_size         = 2,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)


big_plot1

#################################
# phylogenetic diversity
###############################
tree <- read.tree("data/tree.nwk")
tree
comu <- t(comm)
dim(comu)
name.check(tree, comu) #checking names
prunedphy <- prune.sample(t(comu),tree) #cutting tree
prunedphy
plotTree(prunedphy) #visualizing tree

## Converting tree to a distance matrix (cophenetic distance) and performing a PCoA
phydist <- cophenetic(prunedphy)
pcoa <- pcoa(phydist)
barplot(pcoa$values$Eigenvalues)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector <- pcoa$vectors[,1:3]

#convex hull
c <- convhulln(vector,options="Fa", output.options=TRUE)
vol <- c$vol
vol #phylogenetic diversity/volume

f <- hull.build(comm, vector,axes=0)
cont <- hull.contribution(f, relative=TRUE) 
cont #Contribution of each species or individual to the total volume of one or more convex hulls.

##############################################
##### plotting variation in diversity indices
##############################################

data <- read.table("./data/data.txt", header = T)
head(data)
#data with diversities values for each ecoregion according different
#extinction levels. 

#functional
data %>%
  filter(scenario %in% "85_land_61") %>%
  ggplot() +
  aes(x = extinction, y = int_space, colour = ecoregion) +
  geom_line(size = 1.5)+
  geom_point(size = 1.8)+
  theme_classic(base_size = 18) +
  #ylim(0, 0.01)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Interaction diversity",
  ) + 
  guides(fill=guide_legend(title="Ecoregion"))+
  scale_color_manual(name="Ecoregion",
                     labels=c("Atlantic Forest","Cerrado","Chaco", "Pantanal"),values=wes_palette(n=4, name="GrandBudapest2"))


#functional
data %>%
  filter(scenario %in% "85_land_61") %>%
  ggplot() +
  aes(x = extinction, y = fun_space, colour = ecoregion) +
  geom_line(size = 1.5)+
  geom_point(size = 1.8)+
  theme_classic(base_size = 18) +
  #ylim(0, 0.01)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Functional diversity",
  ) + 
  guides(fill=guide_legend(title="Ecoregion"))+
  scale_color_manual(name="Ecoregion",
                     labels=c("Atlantic Forest","Cerrado","Chaco", "Pantanal"),values=wes_palette(n=4, name="GrandBudapest2"))


#phylogenetic
data %>%
  filter(scenario %in% "85_land_61") %>%
  ggplot() +
  aes(x = extinction, y = phy_space, colour = ecoregion) +
  geom_line(size = 1.5)+
  geom_point(size = 1.8)+
  theme_classic(base_size = 18) +
  ylim(0, 0.01)+
  scale_x_continuous(breaks = seq(0, 75, by = 25))+
  labs(
    x = "Number of extinctions (%)",
    y = "Phylogenetic diversity",
  ) + 
  guides(fill=guide_legend(title="Ecoregion"))+
  scale_color_manual(name="Ecoregion",
                     labels=c("Atlantic Forest","Cerrado","Chaco", "Pantanal"),values=wes_palette(n=4, name="GrandBudapest2"))


######################
### network structure
######################

#connectance
con <- networklevel(int, index="connectance")
con

#weighted NODF
nodf.int <- networklevel(int, index="weighted NODF")
nodf.int #nodf
int.random <- nullmodel(int,N=1000, method="vaznull")
nodf.ram.int <- unlist(sapply(int.random, networklevel, index="weighted NODF"))
nodf.ram.int.mean <- mean(nodf.ram.int) 
nodf.ram.int.sd <- sd(nodf.ram.int) 
delta.nodf.int <- nodf.int-nodf.ram.int.mean
zscore.nodf <- delta.nodf.int/nodf.ram.int.sd
zscore.nodf #zscore


#modularity
mod <- computeModules(int, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                    steps = 1E7, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
mod@likelihood #modularity
int.random <- nullmodel(int,N=100, method="vaznull")
mod.ram.int<- sapply(int.random, computeModules, method="Beckett")
mod.ram <- sapply(mod.ram.int, function(x) x@likelihood)
mod.ram.int.mean <- mean(mod.ram) 
mod.ram.int.sd <- sd(mod.ram) 
delta.mod.int <- mod@likelihood-mod.ram.int.mean
zscore.mod <- delta.mod.int/mod.ram.int.sd
zscore.mod #zscore

#H2
h2.int <- networklevel(int, index="H2")
h2.int #specialization
int.random <- nullmodel(int,N=1000, method="vaznull")
h2.ram.int <- unlist(sapply(int.random, networklevel, index="H2"))
h2.ram.int.mean <- mean(h2.ram.int) 
h2.ram.int.sd <- sd(h2.ram.int) 
delta.h2.int <- (h2.int-h2.ram.int.mean)
zscore.h2 <- delta.h2.int/h2.ram.int.sd
zscore.h2 #zscore

#functional complementarity
fc.net <- networklevel(int, index="functional complementarity", dist ='horn', weighted=TRUE)
fc.net

#weighted degree

ed.list <- web2edges(int, both.directions = TRUE,  return = TRUE)
deg <- degree_w(ed.list, alpha=0.5) 
deg <- as.data.frame(deg[1:nrow(int),])
deg

###########################################
# plotting projected range shift
#########################################

data <- read.table("./data/loss.txt", header=T)
data

pant <- data %>%
  filter(biome %in% "pantanal") %>%
  arrange(loss) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(species=factor(species, levels=species)) %>%
  ggplot() +
  aes(x = species, fill = degree, weight = sort(loss)) +
  geom_bar(alpha=1) +
  ylim(-100, +50)+
  labs(
    x = "Species",
    y = "Range shift (%)" ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=90,hjust=1)) 

cerrado <- data %>%
  filter(biome %in% "cerrado") %>%
  arrange(loss) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(species=factor(species, levels=species)) %>%
  ggplot() +
  aes(x = species, fill = degree, weight = sort(loss)) +
  geom_bar(alpha=1) +
  ylim(-100, +50)+
  labs(
    x = "Species",
    y = "Range shift (%)" ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=90,hjust=1)) 

af <- data %>%
  filter(biome %in% "atlantic_forest") %>%
  arrange(loss) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(species=factor(species, levels=species)) %>%
  ggplot() +
  aes(x = species, fill = degree, weight = sort(loss)) +
  geom_bar(alpha=1) +
  ylim(-100, +50)+
  labs(
    x = "Species",
    y = "Range shift (%)" ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=90,hjust=1)) 

chaco <- data %>%
  filter(biome %in% "chaco") %>%
  arrange(loss) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(species=factor(species, levels=species)) %>%
  ggplot() +
  aes(x = species, fill = degree, weight = sort(loss)) +
  geom_bar(alpha=1) +
  ylim(-100, +50)+
  labs(
    x = "Species",
    y = "Range shift (%)" ) +
  scale_fill_distiller(palette = "Greens", direction = 1) +
  theme_classic(base_size = 18)+
  theme(axis.text.x=element_text(angle=90,hjust=1)) 


plot_grid(af, chaco, cerrado, pant, labels=c("A)", "B)", "C)", "D"))


