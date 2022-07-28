# calculate population and prevalence for each country -------------------------
# set-up
library(sf)
library(tidyverse)
library(fuzzyjoin)
library(rgdal)
library(terra)
library(raster)
library(spatial)
library(exactextractr)

# Large files need to be downloaded and saved locally first

datadir <- # set to your data directory #
  
# Files to download and save in datadir:  

# For population counts, use WorldPop's top-down unconstrained population mapping
# (https://www.worldpop.org/methods/top_down_constrained_vs_unconstrained) 
# as obtained through the unconstrained global mosaic 2020 (1km resolution)
# (https://www.worldpop.org/geodata/summary?id=24777).
# 
# Malaria data come from the malaria atlas project
# (https://malariaatlas.org/explorer/#/) raster files called 
# "*Plasmodium falciparum* PR2-10 version 2020", "*Plasmodium vivax* PR1-99 version 2020, 
# "Plasmodium vivax Spatial Limits 2010", and "*Plasmodium falciparum* Spatial Limits 2010."
# 
# Administrative boundaries for level 0 come from https://www.geoboundaries.org/downloadCGAZ.html:
# download shapefile and convert to RDS file using 'sf'

# read in administrative boundaries
admin0 <- readRDS(file.path(datadir, "admin0.RDS"))

# read in Pf and Pv transmission spatial limits and worldpop files 2000 and 2019/2020
spatial_limit_f <- rast(file.path(datadir, 'spatial data', '2010_Pf_Limits_Decompressed.geotiff'))
spatial_limit_v <- rast(file.path(datadir, 'spatial data', '2010_Pv_limits_5k_Decompressed.geotiff'))
worldpop2020 <- rast(file.path(datadir, 'spatial data', 'ppp_2020_1km_Aggregated.tif'))
worldpop2000 <- rast(file.path(datadir, 'spatial data', 'ppp_2000_1km_Aggregated.tif'))

crs(spatial_limit_v); crs(spatial_limit_v); crs(worldpop2000); ext(spatial_limit_f); ext(spatial_limit_v); ext(worldpop2000) # check that files have the same coordinate system: WGS84, EPSG 4326

e <- extent(-180, 180, -60, 84) # set an extent to use for files

# crop spatial limit and worldpop datasets to match so that we can use mask() later on
ext(worldpop2020) <- round(ext(worldpop2020))
ext(worldpop2000) <- round(ext(worldpop2000))

spatial_limit_f <- crop(spatial_limit_f, e)
spatial_limit_v <- crop(spatial_limit_v, e)
worldpop2020 <- crop(worldpop2020, e)
worldpop2000 <- crop(worldpop2000, e)

# the worldpop 2000 dataset uses a large number to indicate missing values
# change these values to missing
worldpop2000[worldpop2000 > 3e+38] = NA
summary(worldpop2000)

# mask population values which are outside of the Pf and Pv transmission areas
malaria_areas2020_pf <- terra::mask(worldpop2020, spatial_limit_f, inverse = F, maskvalues = NA, updatevalue = NA)
malaria_areas2000_pf <- terra::mask(worldpop2000, spatial_limit_f, inverse = F, maskvalues = NA, updatevalue = NA)
malaria_areas2020_pv <- terra::mask(worldpop2020, spatial_limit_v, inverse = F, maskvalues = NA, updatevalue = NA)
malaria_areas2000_pv <- terra::mask(worldpop2000, spatial_limit_v, inverse = F, maskvalues = NA, updatevalue = NA)

# quick visualization
par(mfrow=c(2,1))
raster::plot(malaria_areas2000_pf, breaks = c(0,100,500,1000), main = 'Pf spatial limits + worldpop')
raster::plot(malaria_areas2000_pv, breaks = c(0,100,500,1000), main = 'Pv spatial limits + worldpop')

# write raster and read back in to get file into raster format
makerast <- function(data){
  df.name <- deparse(substitute(data))
  writeRaster(data, file.path(datadir, paste0(df.name, ".tiff")), overwrite=TRUE)
}

makerast(malaria_areas2000_pf)
makerast(malaria_areas2020_pf)
makerast(malaria_areas2000_pv)
makerast(malaria_areas2020_pv)

malaria_areas2000_pf <- rast(file.path(datadir, "malaria_areas2000_pf.tiff"))
malaria_areas2020_pf <- rast(file.path(datadir, "malaria_areas2020_pf.tiff"))
malaria_areas2000_pv <- rast(file.path(datadir, "malaria_areas2000_pv.tiff"))
malaria_areas2020_pv <- rast(file.path(datadir, "malaria_areas2020_pv.tiff"))

# extract population and prevalence values from admin areas within Pf and Pv transmission areas
extractpop <- function(data, data2, species, year){
  # extract population count 
  pop_admins0 <- exact_extract(raster::raster(data), admin0, fun='sum')
  
  # extract PR 
  PR_admins0 <- exact_extract(data2, admin0, fun='mean')

  # regions & countries & ISO names
  ISOregion <- read_csv(file.path(datadir, "ISOregion.csv"))
  
  # bind results into one dataframe
  admins0_all <- admin0 %>% # country boundaries
                 cbind(pop_admins0) %>%  # population
                 cbind(PR_admins0) %>% # PfPR2-10
    # turning 0s to NAs and modifying Philippines name error in dataset
    mutate(PR_admins0 = case_when(PR_admins0==0~NA_real_, TRUE~PR_admins0),
           ISO = case_when(Country=="PHL"~"PHL", TRUE~ISO),
           Country = case_when(Country=="PHL"~"Philippines", TRUE ~ Country)) %>% 
    left_join(ISOregion, by="ISO")
  
  saveRDS(admins0_all, file.path(datadir, paste0("admins0_", year, species, ".rds")))
}

# import PR rasters
PfPR2000 <- raster::raster(file.path(datadir,"/spatial data/2020_GBD2019_Global_PfPR_2000.tif"))
PvPR2000 <- raster::raster(file.path(datadir,"/spatial data/2020_GBD2019_Global_PvPR_2000.tif"))

# run function
extractpop(malaria_areas2000_pf, PfPR2000, '_pf', 2000)
extractpop(malaria_areas2000_pv, PvPR2000, '_pv', 2000)

# read-in a few results and visualize
admins0_2000_pv <- readRDS(file.path(datadir,"admins0_2000_pv.rds"))
admins0_2000_pf <- readRDS(file.path(datadir,"admins0_2000_pf.rds"))
str(admins0_2000_pf)

par(mfrow=c(2,1))
plot(admins0_2000_pf["PR_admins0"], border = 'grey', axes = F, main = "P. falciparum", breaks = seq(0,0.8,0.1))
plot(admins0_2000_pv["PR_admins0"], border = 'grey', axes = F, main = "P. vivax", breaks = seq(0,0.8,0.1))

# there is some bleed-over of PR into neighboring admin0s, perhaps from slightly mismatching shapefiles. E.g. vivax in Kenya. remove cases from countries neighboring endemic countries

# read in MAP excel PR sheets
MAPpf <- read.csv(file.path(datadir, "spatial data/00_PfPR_table_Global_admin0_2000-2019.csv")) %>% filter(Year==2000) %>% dplyr::select(ISO, PfPR_rmean)

MAPpv <- read.csv(file.path(datadir, "spatial data/00_PvPR_table_Global_admin0_2000-2019.csv")) %>% filter(Year==2000) %>% dplyr::select(ISO, PvPR_rmean)

# filter out countries with PfPR and/or Pvpr == 0
admin0_pf <- admins0_2000_pf %>% 
  st_drop_geometry() %>% 
  left_join(MAPpf, by='ISO') %>%
  filter(!(PR_admins0==0 | PfPR_rmean==0 | is.na(PR_admins0) | is.na(PfPR_rmean))) %>%
  dplyr::select(-PfPR_rmean) %>%
  rename(PR_admins0_pf=PR_admins0,
         pop_pf = pop_admins0)

admin0_pv <- admins0_2000_pv %>% 
  st_drop_geometry() %>% 
  left_join(MAPpv, by='ISO') %>%
  filter(!(PR_admins0==0 | PvPR_rmean==0 | is.na(PR_admins0) | is.na(PvPR_rmean))) %>%
  dplyr::select(-PvPR_rmean) %>%
  rename(PR_admins0_pv=PR_admins0,
         pop_pv = pop_admins0)

# save
saveRDS(admin0_pf, file.path(datadir,"admins0_all_2000_pf.rds"))
saveRDS(admin0_pv, file.path(datadir,"admins0_all_2000_pv.rds"))

