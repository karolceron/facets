#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires


#####################################
### land use mask #############
######################################

# Loading packages
library(raster)
library(rgdal)
library(reshape2)
library(letsR)


sp.names <- read.table('sp_names.csv', header = T, sep = ";")
sp.name <- unique(sp.names$species)[10] #choose the species of interest
abr <- unique(sp.names$abr)[10]

# CLIPPING EMNS BY LAND-USE DATA
#### anuran models
refn <- list.files(pattern = ".img")
ref <- raster(refn[1]) # Template for resolution
ref@crs <- CRS("+proj=longlat +datum=WGS84")
plot(ref)

# Land-use data (Li et al. 2017) 
l50A1B <- raster("./data/landuse/world_A1B_2050.tif") #Mitigation
l50A1B <- projectRaster(l50A1B, ref, method="ngb")
plot(l50A1B)

l50A2 <- raster("./data/landuse/world_A2_2050.tif") #BAU
l50A2 <- projectRaster(l50A2, ref, method="ngb")
plot(l50A2)

# Habitat information for anuran species from IUCN
habitat <- read.csv("./data/landuse/tabela_habitat.csv")
head(habitat)

sp_code <- read.csv("./data/landuse/tabela_nomes_sp.csv")
head(sp_code)
sp_h <- merge(sp_code, habitat)

sp_h1 <- data.frame(sp_h$abr, sp_h$X, sp_h$X.1, sp_h$X.2, sp_h$X.2, sp_h$X.4, sp_h$X.5)
sp_h1 <- melt(sp_h1, id = "sp_h.abr")
head(sp_h1)
sp_h1 <- sp_h1[order(sp_h1$sp_h.abr), ]
sp_h1 <- data.frame(sp_h1$sp_h.abr, sp_h1$value)
colnames(sp_h1) <- c("sp_code", "habitat")

dup <- duplicated(sp_h1)
sp_h2 <- sp_h1[-which(dup==TRUE), ]
sp_h2 <- sp_h2[-which(sp_h2$habitat == ""), ]

unique(sp_h2$habitat)

# converting IUCN habitat codes into Xia Li (2017)'s land use classification
# 1 - water
# 2 - forest
# 3 - grassland
# 4 - farmland
# 5 - urban
# 6 - barren

#sp_h2[sp_h2$habitat=="Forest", "habitat"] <- "forest"
sp_h2[sp_h2$habitat==" Wetlands (inland)", "habitat"] <- "water"  
sp_h2[sp_h2$habitat==" Shrubland", "habitat"] <- "grassland"
sp_h2[sp_h2$habitat==" Grassland", "habitat"] <- "grassland"
sp_h2[sp_h2$habitat==" Savanna", "habitat"] <- "grassland"
sp_h2[sp_h2$habitat==" Artificial/Aquatic & Marine", "habitat"] <- "water"
sp_h2[sp_h2$habitat==" Artificial/Terrestrial", "habitat"] <- "farmland"

#write.csv(unique(sp_h2$sp))
#write.csv(sp_h2)

sp_h2[sp_h2$habitat=="forest", "habitat"] <- 2
sp_h2[sp_h2$habitat=="grassland", "habitat"] <- 3
sp_h2[sp_h2$habitat=="farmland", "habitat"] <- 4
sp_h2[sp_h2$habitat=="urban", "habitat"] <- 5
sp_h2[sp_h2$habitat=="barren", "habitat"] <- 6
sp_h2[sp_h2$habitat=="water", "habitat"] <- 1

sp_h2 <- sp_h2[order(sp_h2$sp), ] #species per habitat

sp <- unique(sp_h2$sp)

# ----
# Remove non-habitat
#@45
a <- list.files(pattern=".img")
a
a<-a[c(1,5)] #only img images from 4.5 scenario

`%!in%` = Negate(`%in%`)

     # Ecological niche models 
    
    a1 <- grep(as.character(abr), a)
    a2 <- stack(a[a1])
    names(a2) = a
    plot(a2)
   
    l1 <- crop(l50A1B, a2)             
    plot(l1)
    
    h1 <- sp_h2[sp_h2$sp==abr, "habitat"]
    h1 <- unique(as.numeric(h1))
    
    l1[l1 %!in% h1] <- 0
    l1[l1 %in% h1] <- 1
    
    plot(l1)
    
    l2 <- l1 + a2
    plot(l2)
    
    l2[l2 < 2] <- 0
    l2[l2 == 2] <- 1
    
    names(l2) <- names(a2)
    plot(l2)
    writeRaster(l2, filename=names(l2), bylayer=TRUE,format="GTiff", overwrite=T)



#@85
a <- list.files(pattern=".img")
a
a<-a[c(3,7)] #only img images from 8.5 scenario
a

`%!in%` = Negate(`%in%`)

    a1 <- grep(as.character(abr), a)
    a2 <- stack(a[a1])
    names(a2) = a
    plot(a2)
    
    l1 <- crop(l50A2, a2)             
    plot(l1)
    
    h1 <- sp_h2[sp_h2$sp==sp[38], "habitat"]
    h1 <- unique(as.numeric(h1))
    
    l1[l1 %!in% h1] <- 0
    l1[l1 %in% h1] <- 1
    
    plot(l1)
    
    l2 <- l1 + a2
    plot(l2)
    
    l2[l2 < 2] <- 0
    l2[l2 == 2] <- 1
    
    names(l2) <- names(a2)
    plot(l2)
    writeRaster(l2, filename=names(l2), bylayer=TRUE,format="GTiff", overwrite=T)
  

d <- list.files(pattern =".tif")
d
d<-d[c(1,3)] #only .tif images of 4.5 scenario
d
d1 <- stack(d)
names(d1) <- d
d2 <- as.data.frame(d1)

head(d2)
d3 <- colSums(d2, na.rm = T)
d4 <- data.frame(d3)
d4$scenario <- "4.5"
colnames(d4) <- c("cells_after_mask", "scenario")
d4
write.csv(d4)


e <- list.files(pattern =".tif")
e<-e[c(1,3)] #only .tif images
e1 <- stack(e)
names(e1) <- e
e2 <- as.data.frame(e1)

head(e2)
e3 <- colSums(e2, na.rm = T)
e4 <- data.frame(e3)
e4$scenario <- "8.5"

colnames(e4) <- c("cells_after_mask", "scenario")
e4
write.csv(e4)


