#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes in predator-prey networks"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires


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

#read data
int <- t(read.table("data/pant.txt"))
int<-as.data.frame(int)
rowSums(int) #checking data
colSums(int) #checking data

int<-int[, which(colSums(int) != 0)] #removing zero columns
int


int$Formicidae<-log(1+int$Formicidae) #transforming more abundant data
int$Isoptera<-log(1+int$Isoptera) #transforming more abundant data

######################################
# interaction diversity 
########################################

#creating community matrix
anura<-row.names(int)
site<-rep(1, nrow(int)) 
comm <- data.frame(matrix(nrow = 1,data = site))
colnames(comm) <- anura
rownames(comm) <- 'site'
comm

## Converting interaction to a distance matrix (horn distance) and transforming them with a PCoA
horn=vegdist(int, method="horn")
pcoa=pcoa(horn)
barplot(pcoa$values$Eigenvalues)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector=pcoa$vectors[,1:3]

#convex hull 
c=convhulln(vector,options="Fa", output.options=TRUE)
vol=c$vol #volume
vol #interaction diversity/volume
f<-hull.build(comm, vector,axes=0)
cont<-hull.contribution(f, relative=TRUE) 
cont #Contribution of each species or individual to the total volume of one or 
#more convex hulls.


######################################
# functional diversity 
########################################

traits <- read.table("data/traits.txt", header = T) #read data
traits<-as.data.frame(traits)

traits <- traits[which(rownames(traits) %in% colnames(comm)), ] #cutting traits according community 
traits
dim(traits)

traits<-traits[, which(colSums(traits) != 0)] #removing zero columns
traits

## Converting traits to a distance matrix (gower distance) and transforming them with a PCoA
hab=prep.binary(traits[,4:6],col.blocks=3,label="habitat")
mass=data.frame(body=traits$mass)
clutch=data.frame(clutch=traits$clutch_size)
size=data.frame(size=traits$size)
mode=data.frame(mode=traits$mode)
ktab=ktab.list.df(list(size, clutch, mass,hab, mode), rownames = rownames(traits))
dist=dist.ktab(ktab, type= c("Q", "Q", "Q","B", "O")) #gower mixed
pcoa=pcoa(dist)
barplot(pcoa$values$Eigenvalues)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector=pcoa$vectors[,1:3]

#convex hull 
c=convhulln(vector,options="Fa", output.options=TRUE)
vol=c$vol
vol #functional diversity/volume
f<-hull.build(comm, vector,axes=0)
cont<-hull.contribution(f, relative=TRUE) 
cont #Contribution of each species or individual to the total volume of one or more convex hulls.

#################################
# phylogenetic diversity
###############################
tree=read.tree("data/tree.nwk")
tree
comu<-t(comm)
dim(comu)
name.check(tree, comu) #checking names
prunedphy <-prune.sample(t(comu),tree) #cutting tree
prunedphy
plotTree(prunedphy) #visualizing tree

## Converting tree to a distance matrix (cophenetic distance) and transforming them with a PCoA
phydist=cophenetic(prunedphy)
pcoa=pcoa(phydist)
barplot(pcoa$values$Eigenvalues)
(pcoa$values$Eigenvalues[1]+pcoa$values$Eigenvalues[2]+pcoa$values$Eigenvalues[3])/sum(pcoa$values$Eigenvalues[pcoa$values$Eigenvalues> 0])

#saving vectors
vector=pcoa$vectors[,1:3]

#convex hull
c=convhulln(vector,options="Fa", output.options=TRUE)
vol=c$vol
vol #phylogenetic diversity/volume

f<-hull.build(comm, vector,axes=0)
cont<-hull.contribution(f, relative=TRUE) 
cont #Contribution of each species or individual to the total volume of one or more convex hulls.


######################
### network structure
######################

#connectance
con.suarede <-networklevel(suarede, index="connectance")
con.suarede

#weighted NODF
nodf.suarede <-networklevel(suarede, index="weighted NODF")
nodf.suarede #nodf
suarede.random <- nullmodel(suarede,N=1000, method="vaznull")
nodf.ram.suarede <- unlist(sapply(suarede.random, networklevel, index="weighted NODF"))
nodf.ram.suarede.mean <- mean(nodf.ram.suarede) 
nodf.ram.suarede.sd <- sd(nodf.ram.suarede) 
delta.nodf.suarede <- nodf.suarede-nodf.ram.suarede.mean
zscore.nodf<-delta.nodf.suarede/nodf.ram.suarede.sd
zscore.nodf #zscore


#modularity
modsuarede<-computeModules(suarede, method="Beckett", deep = FALSE, deleteOriginalFiles = TRUE, 
                           steps = 1E7, tolerance = 1e-10, experimental = FALSE, forceLPA=FALSE)
modsuarede@likelihood #modularity
suarede.random <- nullmodel(suarede,N=100, method="vaznull")
mod.ram.suarede<- sapply(suarede.random, computeModules, method="Beckett")
mod.ram <- sapply(mod.ram.suarede, function(x) x@likelihood)
mod.ram.suarede.mean <- mean(mod.ram) 
mod.ram.suarede.sd <- sd(mod.ram) 
delta.mod.suarede <- modsuarede@likelihood-mod.ram.suarede.mean
zscore.mod<-delta.mod.suarede/mod.ram.suarede.sd
zscore.mod #zscore

#H2
h2.suarede <-networklevel(suarede, index="H2")
h2.suarede #specialization
suarede.random <- nullmodel(suarede,N=1000, method="r2dtable")
h2.ram.suarede <- unlist(sapply(suarede.random, networklevel, index="H2"))
h2.ram.suarede.mean <- mean(h2.ram.suarede) 
h2.ram.suarede.sd <- sd(h2.ram.suarede) 
delta.h2.suarede <- (h2.suarede-h2.ram.suarede.mean)
zscore.h2<-delta.h2.suarede/h2.ram.suarede.sd
zscore.h2 #zscore

#functional complementarity
fc.suarede <-networklevel(suarede, index="functional complementarity", dist ='horn', weighted=TRUE)
fc.suarede

#weighted degree
deg<-degree_w(suarede,measure=c("degree"), alpha=0.5) 
deg<-as.data.frame(deg)
deg


