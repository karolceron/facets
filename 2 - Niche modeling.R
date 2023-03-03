#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires


## Warning ##

# After downloading occurrence data from script #1, data was manually curated
# to avoid wrong distribution points and to include data occurrences from
# the points sampled

# Curated data are available in folder" occurrence_data"



#####################################
### niche modeling #############
######################################

# please note that the niche modeling analysis was performed using biomod2 version 4.1
# as biomod2 package is under constantly modification, as biomod2 package is 
# under constant modification, performing these analyses using the
# updated version may do not work.

#to run the maxent model please see instructions
# (https://r-forge.r-project.org/scm/viewvc.php/*checkout*/pkg/biomod2/inst/doc/Include_MAXENT.pdf?revision=367&root=biomod) for more information



## install the 4.1-2 version of biomod2
install.packages("./package/biomod2-4.1-1.tar.gz", repo=NULL, type="source")

#If you have trouble with this part refer to:
#https://biomodhub.github.io/biomod2/articles/examples_1_mainFunctions.html

## load the required packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(sf)
library(rgdal)
library(ncdf4)
library(spThin)
library(usdm)

## read data ----
setwd("D:/Dropbox/Backup/Orientacao/Karol_Ceron/Mathias-Karol/codes/MS-foodwebs-EcolLett/code")

sp.names <- read.table('sp_names.csv', header = T, sep = ";")
sp.name <- unique(sp.names$species)[10] #choose the species of interest
abr <- unique(sp.names$abr)[10] #indicate the abbreviation of species 


ProLau_occ <- read.table(paste0("./occurrence_data/",abr, '.txt'), sep=",", header = T)
summary(ProLau_occ)

#thinning data points
occ_thin <- thin(ProLau_occ, verbose=T, 
               long.col= "long", 
               lat.col="lat", 
               spec.col= abr,  
               thin.par= 3, 
               reps=1,
               locs.thinned.list.return=T,
               write.files=F)

occ_thin1 <- data.frame(occ_thin)
occ_thin <- data.frame(occ_thin)
occ_thin$sp1 <- 1

# read files from  Bioclim
files <- list.files('./data/bio/', pattern='.tif$', full.names=TRUE)
predictors <- stack(files)

# crop predictors
e <- extent(-90,-30,-50,15)
layers <- crop(predictors,  e)

#plotting layers
plot(layers$elev) #example with elevation
points(occ_thin[,1:2])

# remove variables that are highly correlated
spx <- extract(layers, occ_thin1)
spx <- data.frame(spx)
v <- vifstep(spx, th=5) 
v
bioclim <- exclude(layers, v)
bioclim <- stack(bioclim)


## format the data 
niche_data <- 
  BIOMOD_FormatingData(
    resp.var = occ_thin['sp1'],
    resp.xy = occ_thin[, c('Longitude','Latitude')],
    expl.var = bioclim,
    resp.name = abr,
    PA.nb.rep = 2,
    PA.nb.absences = 500,
    PA.strategy = 'random'
  )

## formatted object summary
niche_data


## define individual models options ---- 
ProLau_opt <- 
  BIOMOD_ModelingOptions(
    GLM = list(type = 'quadratic', interaction.level = 1),
    GBM = list(n.trees = 1000),
    GAM = list(algo = 'GAM_mgcv'),
    CTA = ,
    SRE = ,
    FDA = ,
    MARS = list(interaction.level = 1),
    RF=,
    #MAXENT.Phillips =
  )

## run the individual models ----
niche_models <- 
  BIOMOD_Modeling(
    bm.format = niche_data,
    models = c("GLM", "GBM", "RF", "GAM", "CTA", "SRE", "FDA", "MARS"),#, "MAXENT.Phillips"),
    bm.options = ProLau_opt,
    nb.rep = 2,
    data.split.perc = 80,
    var.import = 1,
    modeling.id = "demo1")


## assess individual models quality ----
## get models evaluation scores
niche_models_scores <- get_evaluations(niche_models)
niche_models_scores


## check variable importance
niche_models_var_import <- get_variables_importance(niche_models)

## compute the mean of variable importance by algorithm
apply(niche_models_var_import, c(1,2), mean, na.rm = TRUE)

## run the ensemble models ----
ensemble_models <- 
  BIOMOD_EnsembleModeling(
    bm.mod = niche_models,
    em.by = 'all', #all models
    metric.select = c('TSS'),
    metric.eval = c('TSS', 'ROC'),
    metric.select.thresh =c(0.8), 
    prob.mean = FALSE,
    prob.cv = TRUE, 
    committee.averaging = TRUE,
    prob.mean.weight = TRUE,
    var.import = 0)


## asses ensemble models quality ----
(ensemble_models_scores <- get_evaluations(ensemble_models))

## do models projections ----
## current projections
niche_models_proj_current <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    proj.name = "current",
    build.clamping.mask = TRUE,
    compress = FALSE)


ensemble_models_proj_current <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_current,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)



#################################################################
## future projections scenario_4.5 
#2041-2060
## load BCC-CSM2-MR variables
bioclim_BCC45_2041 <-
  stack('./data/4.5/2041/BCC-CSM2-MR_ssp245_2041-2060.tif')

names(bioclim_BCC45_2041)

names(bioclim_BCC45_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                               , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                               , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                               , 'bio_19')
names(bioclim_BCC45_2041)

bioclim_BCC45_2041 <- crop(bioclim_BCC45_2041,  e)
bioclim_BCC45_2041 <- stack(bioclim_BCC45_2041, bioclim$elev)


#selected variables

bioclim_BCC45_2041 <- subset(bioclim_BCC45_2041, v@results$Variables) 
bioclim_BCC45_2041 <- stack(bioclim_BCC45_2041)


niche_models_proj_BCC45_2041 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_BCC45_2041,
    proj.name = "BCC45_2041",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)


ensemble_models_proj_BCC45_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_BCC45_2041,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)



## load 2061 bioclim variables

bioclim_BCC45_2061 <-
  stack('./data/4.5/2061/bioc_BCC-CSM2-MR_ssp245_2061-2080.tif')

names(bioclim_BCC45_2061)

names(bioclim_BCC45_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                               , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                               , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                               , 'bio_19')
names(bioclim_BCC45_2061)


bioclim_BCC45_2061 <- crop(bioclim_BCC45_2061,  e)
bioclim_BCC45_2061 <- stack(bioclim_BCC45_2061, bioclim$elev)


#selected variables

bioclim_BCC45_2061<-subset(bioclim_BCC45_2061, v@results$Variables) 
bioclim_BCC45_2061 <- stack(bioclim_BCC45_2061)


niche_models_proj_2061_BCC45 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_BCC45_2061,
    proj.name = "2061_BC45",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_BC45 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_BCC45,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## compute Species Range Change (SRC) ----
## load binary projections
ProLau_bin_proj_current <- 
  stack( 
    c(wm = paste0(abr,"/proj_current/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))


plot(ProLau_bin_proj_current)

ProLau_bin_proj_2041_BC45 <- 
  stack( 
    c(wm = paste0(abr,"/proj_BCC45_2041/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

plot(ProLau_bin_proj_2041_BC45)

ProLau_bin_proj_2061_BC45 <- 
  stack( 
    c(wm = paste0(abr,"/proj_2061_BC45//individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

plot(ProLau_bin_proj_2061_BC45)

## SRC current -> 2041
SRC_current_2041_BC45 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_BC45)

SRC_current_2041_BC45$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_BC45 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_BC45)

SRC_current_2061_BC45$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_BC45$Diff.By.Pixel, 
    SRC_current_2061_BC45$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf("bicolor_CanESM5.pdf", width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1)
)
#dev.off()

######################################################################
## future projections cenario_4.5 
#2041-2060
## load CanESM5 variables
bioclim_CanESM5_2041 <-
  stack('./data/4.5/2041/CanESM5_ssp245_2041-2060.tif')

names(bioclim_CanESM5_2041)

names(bioclim_CanESM5_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                 , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                 , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                 , 'bio_19')
names(bioclim_CanESM5_2041)


bioclim_CanESM5_2041 <- crop(bioclim_CanESM5_2041,  e)
bioclim_CanESM5_2041 <- stack(bioclim_CanESM5_2041, bioclim$elev)


#selected variables
bioclim_CanESM5_2041 <- subset(bioclim_CanESM5_2041, v@results$Variables) 
bioclim_CanESM5_2041 <- stack(bioclim_CanESM5_2041)


niche_models_proj_CanESM5_2041 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_CanESM5_2041,
    proj.name = "CanESM5_2041",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)


ensemble_models_proj_CanESM5_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_CanESM5_2041,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## load 2061 bioclim variables

bioclim_CanESM5_2061 <-
  stack('./data/4.5/2061/CanESM5_ssp245_2061-2080.tif')

names(bioclim_CanESM5_2061)

names(bioclim_CanESM5_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                 , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                 , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                 , 'bio_19')
names(bioclim_CanESM5_2061)



bioclim_CanESM5_2061 <- crop(bioclim_CanESM5_2061,  e)
bioclim_CanESM5_2061 <- stack(bioclim_CanESM5_2061, bioclim$elev)


#selected variables
bioclim_CanESM5_2061 <- subset(bioclim_CanESM5_2061, v@results$Variables) 
bioclim_CanESM5_2061 <- stack(bioclim_CanESM5_2061)


niche_models_proj_2061_CanESM5 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_CanESM5_2061,
    proj.name = "2061_CanESM5",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_CanESM5 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_CanESM5,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## compute Species Range Change (SRC) ----
## load binary projections

ProLau_bin_proj_2041_CanESM5 <- 
  stack( 
    c(wm = paste0(abr,"/proj_CanESM5_2041/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

ProLau_bin_proj_2061_CanESM5 <- 
  stack( 
    c(wm =  paste0(abr,"/proj_2061_CanESM5/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

## SRC current -> 2041
SRC_current_2041_CanESM5 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_CanESM5)

SRC_current_2041_CanESM5$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_CanESM5 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_CanESM5)

SRC_current_2061_CanESM5$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_CanESM5$Diff.By.Pixel, 
    SRC_current_2061_CanESM5$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf(paste0(sp.name,"_CanESM5.pdf"), width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1)
)
#dev.off()

######################################################################
## future projections cenario_4.5 
#2041-2060
## load MIROC6 variables
bioclim_MIROC6_2041 <-
  stack('./data/4.5/2041/MIROC6_ssp245_2041-2060.tif')

names(bioclim_MIROC6_2041)

names(bioclim_MIROC6_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                , 'bio_19')
names(bioclim_MIROC6_2041)

bioclim_MIROC6_2041 <- crop(bioclim_MIROC6_2041,  e)
bioclim_MIROC6_2041 <- stack(bioclim_MIROC6_2041, bioclim$elev)


#selected variables
bioclim_MIROC6_2041 <- subset(bioclim_MIROC6_2041, v@results$Variables) 
bioclim_MIROC6_2041 <- stack(bioclim_MIROC6_2041)


niche_models_proj_MIROC6_2041 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_MIROC6_2041,
    proj.name = "MIROC6_2041",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_MIROC6_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_MIROC6_2041,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## load 2061 bioclim variables

bioclim_MIROC6_2061 <-
  stack('./data/4.5/2061/MIROC6_ssp245_2061-2080.tif')

names(bioclim_MIROC6_2061)

names(bioclim_MIROC6_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                , 'bio_19')
names(bioclim_MIROC6_2061)

bioclim_MIROC6_2061 <- crop(bioclim_MIROC6_2061,  e)
bioclim_MIROC6_2061 <- stack(bioclim_MIROC6_2061, bioclim$elev)


#selected variables
bioclim_MIROC6_2061 <- subset(bioclim_MIROC6_2061, v@results$Variables) 
bioclim_MIROC6_2061 <- stack(bioclim_MIROC6_2061)



niche_models_proj_2061_MIROC6 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_MIROC6_2061,
    proj.name = "2061_MIROC6",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_MIROC6 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_MIROC6,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## compute Species Range Change (SRC) ----
## load binary projections

ProLau_bin_proj_2041_MIROC6 <- 
  stack( 
    c(wm = paste0(abr,"/proj_MIROC6_2041/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

ProLau_bin_proj_2061_MIROC6 <- 
  stack( 
    c(wm = paste0(abr,"/proj_2061_MIROC6/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

## SRC current -> 2041
SRC_current_2041_MIROC6 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_MIROC6
  )

SRC_current_2041_MIROC6$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_MIROC6 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_MIROC6
  )

SRC_current_2061_MIROC6$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_MIROC6$Diff.By.Pixel, 
    SRC_current_2061_MIROC6$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf(paste0(abr,"_MIROC6.pdf"), width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1))

#dev.off()

###################################################################
##################### 8.5 ##########################################
## future projections cenario_8.8 
#2041-2060
## load BCC-CSM2-MR variables

bioclim_BCC85_2041 <-
  stack('./data/8.5/2041/BCC-CSM2-MR_ssp585_2041-2060.tif')

names(bioclim_BCC85_2041)

names(bioclim_BCC85_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                               , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                               , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                               , 'bio_19')
names(bioclim_BCC85_2041)


bioclim_BCC85_2041 <- crop(bioclim_BCC85_2041,  e)
bioclim_BCC85_2041 <- stack(bioclim_BCC85_2041, bioclim$elev)

#selected variables
bioclim_BCC85_2041 <- subset(bioclim_BCC85_2041, v@results$Variables) 
bioclim_BCC85_2041 <- stack(bioclim_BCC85_2041)

niche_models_proj_BCC85_2041 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_BCC85_2041,
    proj.name = "BCC85_2041",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_BCC85_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_BCC85_2041,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)

## load 2061 bioclim variables

bioclim_BCC85_2061 <-
  stack('./data/8.5/2061/BCC-CSM2-MR_ssp585_2061-2080.tif')

names(bioclim_BCC85_2061)

names(bioclim_BCC85_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                               , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                               , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                               , 'bio_19')
names(bioclim_BCC85_2061)


bioclim_BCC85_2061 <- crop(bioclim_BCC85_2061,  e)
bioclim_BCC85_2061 <- stack(bioclim_BCC85_2061, bioclim$elev)

#selected variables
bioclim_BCC85_2061 <- subset(bioclim_BCC85_2061, v@results$Variables) 
bioclim_BCC85_2061 <- stack(bioclim_BCC85_2061)


niche_models_proj_2061_BCC85 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_BCC85_2061,
    proj.name = "2061_BC85",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_BC85 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_BCC85,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## compute Species Range Change (SRC) ----
## load binary projections

ProLau_bin_proj_2041_BC85 <- 
  stack( 
    c(wm = paste0(abr,"/proj_BCC85_2041/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

ProLau_bin_proj_2061_BC85 <- 
  stack( 
    c(wm = paste0(abr,"/proj_2061_BC85/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))


## SRC current -> 2041
SRC_current_2041_BC85 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_BC85
  )

SRC_current_2041_BC85$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_BC85 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_BC85
  )

SRC_current_2061_BC85$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_BC85$Diff.By.Pixel, 
    SRC_current_2061_BC85$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf(paste0(abr,"_BC85.pdf"), width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1)
)

#dev.off()

######################################################################
## future projections cenario_8.5 
#2041-2060
## load CanESM5 variables
bioclim_CanESM5_2041 <-
  stack('./data/8.5/2041/CanESM5_ssp585_2041-2060.tif')

names(bioclim_CanESM5_2041)

names(bioclim_CanESM5_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                 , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                 , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                 , 'bio_19')
names(bioclim_CanESM5_2041)


bioclim_CanESM5_2041 <- crop(bioclim_CanESM5_2041,  e)
bioclim_CanESM5_2041 <- stack(bioclim_CanESM5_2041, bioclim$elev)

#selected variables
bioclim_CanESM5_2041 <- subset(bioclim_CanESM5_2041, v@results$Variables) 
bioclim_CanESM5_2041 <- stack(bioclim_CanESM5_2041)

niche_models_proj_CanESM5_2041_85 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_CanESM5_2041,
    proj.name = "CanESM5_85_2041",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_CanESM5_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_CanESM5_2041_85,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## load 2061 bioclim variables

bioclim_CanESM5_2061 <-
  stack('./data/8.5/2061/CanESM5_ssp585_2061-2080.tif')

names(bioclim_CanESM5_2061)

names(bioclim_CanESM5_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                 , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                 , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                 , 'bio_19')
names(bioclim_CanESM5_2061)


bioclim_CanESM5_2061 <- crop(bioclim_CanESM5_2061,  e)
bioclim_CanESM5_2061 <- stack(bioclim_CanESM5_2061, bioclim$elev)

#selected variables
bioclim_CanESM5_2061 <- subset(bioclim_CanESM5_2061, v@results$Variables) 
bioclim_CanESM5_2061 <- stack(bioclim_CanESM5_2061)

niche_models_proj_2061_CanESM5_85 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_CanESM5_2061,
    proj.name = "2061_CanESM5_85",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_CanESM5_85 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_CanESM5_85,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)


## compute Species Range Change (SRC) ----
## load binary projections


ProLau_bin_proj_2041_CanESM5 <- 
  stack( 
    c(wm = paste0(abr,"/proj_CanESM5_85_2041/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

ProLau_bin_proj_2061_CanESM5 <- 
  stack( 
    c(wm = paste0(abr,"/proj_2061_CanESM5_85/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

## SRC current -> 2041
SRC_current_2041_CanESM5 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_CanESM5
  )

SRC_current_2041_CanESM5$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_CanESM5 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_CanESM5
  )

SRC_current_2061_CanESM5$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_CanESM5$Diff.By.Pixel, 
    SRC_current_2061_CanESM5$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf(paste0(abr,"_CanESM5_85.pdf"), width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1)
)

#dev.off()

######################################################################
## future projections cenario_8.5 
#2041-2060
## load MIROC6 variables
bioclim_MIROC6_2041 <-
  stack('./data/8.5/2041/MIROC6_ssp585_2041-2060.tif')

names(bioclim_MIROC6_2041)

names(bioclim_MIROC6_2041) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                , 'bio_19')
names(bioclim_MIROC6_2041)


bioclim_MIROC6_2041 <- crop(bioclim_MIROC6_2041,  e)
bioclim_MIROC6_2041 <- stack(bioclim_MIROC6_2041, bioclim$elev)

#selected variables
bioclim_MIROC6_2041 <- subset(bioclim_MIROC6_2041, v@results$Variables) 
bioclim_MIROC6_2041 <- stack(bioclim_MIROC6_2041)

niche_models_proj_MIROC6_2041_85 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_MIROC6_2041,
    proj.name = "MIROC6_2041_85",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_MIROC6_2041 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_MIROC6_2041_85,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)

## load 2061 bioclim variables

bioclim_MIROC6_2061 <-
  stack('./data/8.5/2061/MIROC6_ssp585_2061-2080.tif')

names(bioclim_MIROC6_2061)

names(bioclim_MIROC6_2061) <- c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6'
                                , 'bio_7', 'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12'
                                , 'bio_13', 'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18'
                                , 'bio_19')
names(bioclim_MIROC6_2061)


bioclim_MIROC6_2061 <- crop(bioclim_MIROC6_2061,  e)
bioclim_MIROC6_2061 <- stack(bioclim_MIROC6_2061, bioclim$elev)

#selected variables
bioclim_MIROC6_2061 <- subset(bioclim_MIROC6_2061, v@results$Variables) 
bioclim_MIROC6_2061 <- stack(bioclim_MIROC6_2061)

niche_models_proj_2061_MIROC6_85 <- 
  BIOMOD_Projection(
    bm.mod = niche_models,
    new.env = bioclim_MIROC6_2061,
    proj.name = "2061_MIROC6_85",
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    build.clamping.mask = TRUE,
    compress = FALSE)

ensemble_models_proj_2061_MIROC6 <- 
  BIOMOD_EnsembleForecasting(
    bm.em = ensemble_models,
    bm.proj = niche_models_proj_2061_MIROC6_85,
    models.chosen = 'all',
    metric.binary = 'all',
    metric.filter = 'all',
    output.format = ".img",
    do.stack = FALSE,
    compress = FALSE)

## compute Species Range Change (SRC) ----
## load binary projections


ProLau_bin_proj_2041_MIROC6 <- 
  stack( 
    c(wm = paste0(abr,"/proj_MIROC6_2041_85/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

ProLau_bin_proj_2061_MIROC6 <- 
  stack( 
    c(wm = paste0(abr,"/proj_2061_MIROC6_85/individual_projections/", abr, "_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img")
    ))

## SRC current -> 2041
SRC_current_2041_MIROC6 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2041_MIROC6
  )

SRC_current_2041_MIROC6$Compt.By.Models

## SRC current -> 2061
SRC_current_2061_MIROC6 <- 
  BIOMOD_RangeSize(
    ProLau_bin_proj_current,
    ProLau_bin_proj_2061_MIROC6
  )

SRC_current_2061_MIROC6$Compt.By.Models

ProLau_src_map <- 
  stack(
    SRC_current_2041_MIROC6$Diff.By.Pixel, 
    SRC_current_2061_MIROC6$Diff.By.Pixel
  )
names(ProLau_src_map) <- c("cur-2041", "cur-2061")

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

#pdf(paste0(abr,"_MIROC6_85.pdf"), width = 10, height = 8)

rasterVis::levelplot( 
  ProLau_src_map, 
  main = sp.name,
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,1)
)


#dev.off()


##################################
######### merging ensemble maps
#################################

######################
######## 4.5 ########
#####################


#2041

m1 <- sum(ProLau_bin_proj_2041_BC45, ProLau_bin_proj_2041_CanESM5,ProLau_bin_proj_2041_MIROC6)
plot(m1)

mask1 <- m1 == 3
plot(mask1)

writeRaster(mask1, paste0(abr, ".2041.45.img"))

#2061
m2 <- sum(ProLau_bin_proj_2061_BC45, ProLau_bin_proj_2061_CanESM5,ProLau_bin_proj_2061_MIROC6)
plot(m2)

mask2 <- m2 == 3
plot(mask2)

writeRaster(mask2, paste0(abr, ".2061.45.img"))

######################
######## 8.5 ########
#####################

#2041

m3 <- sum(ProLau_bin_proj_2041_BC85, ProLau_bin_proj_2041_CanESM5,ProLau_bin_proj_2041_MIROC6)
plot(m3)

mask3 <- m3 == 3
plot(mask3)

writeRaster(mask3, paste0(abr, ".2041.85.img"))

# 2061
m4 <- sum(ProLau_bin_proj_2061_BC85, ProLau_bin_proj_2061_CanESM5,ProLau_bin_proj_2061_MIROC6)
plot(m4)

mask4 <- m4 == 3
plot(mask4)

writeRaster(mask4, paste0(abr, ".2061.85.img"))



## compute the number of cells lost in climate scenario

a<-(getValues(ProLau_bin_proj_current))
colSums(a, na.rm=T)


q <- 
  stack( 
    c(wm = paste0(abr, ".2041.45.img")))
b<-(getValues(q))
colSums(b, na.rm=T) 

s <- 
  stack( 
    c(wm = paste0(abr, ".2061.45.img")))
c<-(getValues(s))
colSums(c, na.rm=T)

#### 8.5

qq <- 
  stack( 
    c(wm = paste0(abr, ".2041.85.img")))
b<-(getValues(qq))
colSums(b, na.rm=T)

ss <- 
  stack( 
    c(wm = paste0(abr, ".2061.85.img")))
c<-getValues(ss)
colSums(c, na.rm=T)




