#Example of code used in "Decoupled responses of biodiversity facets driven 
#from anuran vulnerability to climate and land use changes"
# by Karoline Ceron, Lilian P. Sales, Diego J. Santana, and Mathias M. Pires


######## installing and loading packages ############
library(devtools)
#install_github("ropensci/CoordinateCleaner")
library(CoordinateCleaner)
library(countrycode)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)
library(rnaturalearthdata)
#devtools::install_github("liibre/Rocc")
library(Rocc)
library(raster)
library(ggplot2)

#################################
#obtain data from GBIF via rgbif
#################################
sp.names <- read.table('sp_names.csv', header = T, sep = ";")
sp.name <- unique(sp.names$species)[10] #choose the species of interest
abr <- unique(sp.names$abr)[10] #indicate the abbreviation of species 

dat <- occ_search(scientificName = sp.name, hasCoordinate = T)
dat <- dat$data
dim(dat)


#select columns of interest
dat <- dat %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)


# remove records without coordinates
dat <- dat%>%
  filter(!is.na(dat$decimalLongitude))%>%
  filter(!is.na(dat$decimalLatitude))

#plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.5)+
  theme_bw()

#convert country code from ISO2c to ISO3c
dat$countryCode <- countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

#flag problems
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries")) # most test are on by default
summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")


#Exclude problematic records
dat_cl <- dat[flags$.summary,]

#The flagged records
dat_fl <- dat[!flags$.summary,]

#to avoid specifying it in each function
names(dat)[2:3] <- c("decimallongitude", "decimallatitude")

clean <- dat%>%
  cc_val()%>%
  cc_equ()%>%
  cc_cap()%>%
  cc_cen()%>%
  cc_coun(iso3 = "countryCode")%>%
  cc_gbif()%>%
  cc_inst()%>%
  cc_sea()%>%
  cc_zero()%>%
  cc_outl()%>%
  cc_dupl()

dat %>%
  as_tibble() %>% 
  mutate(val = cc_val(., value = "flagged"),
         sea = cc_sea(., value = "flagged"))


#Testing temporal outliers on taxon level
flags <- cf_age(x = dat_cl,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                taxon = "species", 
                min_age = "year", 
                max_age = "year", 
                value = "flagged")

#dat_cl=dat
dat_cl[!flags, "year"]
dat_cl <- dat_cl[flags, ]

#Remove records with low coordinate precision
hist(dat_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)

dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))


#Maintain records with according basis of record
table(dat$basisOfRecord)

dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION" | 
                   basisOfRecord == "OBSERVATION" |
                   basisOfRecord == "PRESERVED_SPECIMEN")

#Individual count
table(dat_cl$individualCount)

dat_cl <- dat_cl%>%
  filter(individualCount > 0 | is.na(individualCount))%>%
  filter(individualCount < 99 | is.na(individualCount)) 

#Age of records
table(dat_cl$year)

dat_cl <- dat_cl%>%
  filter(year > 1980) #removing data before XXXX years

#checking taxonomy
table(dat_cl$family)

#dat_cl <- dat_cl%>%
 # filter(family == 'Microhylidae') ## restrict data to a determined taxonomic family

dim(dat_cl)


#plotting
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat_cl, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.5)+
  theme_bw()

dim(dat_cl)



####################################
##obtain data from specieslink
####################################

splink <- rspeciesLink(filename = "ex01", species =  sp.name,
                       Scope = "animals", Coordinates = "YES")

dim(splink)

#select columns of interest
sp <- splink %>%
  dplyr::select(scientificName, decimalLongitude, decimalLatitude, individualCount,
                family, year, basisOfRecord, institutionCode, country)

# remove records without coordinates
sp <- sp%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

sp <- data.frame(sp)
head(sp)

#to avoid specifying it in each function
names(sp)[2:3] <- c("decimallongitude", "decimallatitude")

#formatting data
sp <- data.frame(sp)
sp1 <- sp[,1:3]
gbif <- dat_cl[,1:3]

names(sp1)[1:3] <- c(abr, "Longitude", "Latitude")
names(gbif) <- names(sp1) 

#merging gbid and specislink data
final <- rbind(gbif, sp1)
dim(final)
final

distinct_data <- dplyr::distinct(final)
dim(distinct_data)

distinct_data[,abr] <- 1


#plotting
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = distinct_data, aes(x = as.numeric(distinct_data$Longitude), y = as.numeric(distinct_data$Latitude)),
             colour = "darkred", size = 0.5)+
  theme_bw()

#saving data on a .csv file
write.csv(distinct_data, paste0(abr, '.csv'), row.names = FALSE)

