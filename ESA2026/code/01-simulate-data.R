
# This tutorial uses SpFut.flexiSDMv1.1.2. Using it with
# a different version may cause errors. For instructions on how
# to install SpFut.flexiSDMv1.1.2, see 00-set-up.R



# This script simulates the Great Salt Lake Cyclops distribution
# and four datasets.


# Set up environment ----
library(SpFut.flexiSDM)
library(SpFut.covariates)
library(sf)
library(tidyverse)
library(rnaturalearth)

# rnaturalearth requires rnaturalearthhires, which is not on github. If
# you do not already have it, install it
remotes::install_github("ropensci/rnaturalearthhires")


# set seed so results are reproducible
set.seed(1)





# Region ----


# Get some useful landmarks to define the region

# Great Salt Lake
GSL <- ne_download(scale = 10, type = "lakes",
                   category = "physical", returnclass = "sf") %>%
  filter(name == "Great Salt Lake")

# Salt Lake City point
slc_point <- ne_download(scale = 10, type = "populated_places",
                         category = "cultural", returnclass = "sf") %>%
  filter(NAME == "Salt Lake City")

# Salt Lake City outline: the urban area polygon that contains the SLC point
urban <- ne_download(scale = 10, type = "urban_areas",
                     category = "cultural", returnclass = "sf") %>%
  st_transform(st_crs(slc_point))
SLC <- urban[st_intersects(urban, slc_point, sparse = FALSE), ]

ggplot() +
  geom_sf(data = GSL, fill = "blue") +
  geom_sf(data = SLC, fill = "brown")



# use Great Salt Lake as starting area for region, add buffer
buffer <- 100000
gridstart <- GSL %>% 
  st_transform(crs = 3857) %>% 
  st_buffer(buffer)

# create grid
cellarea <- 25e6  # hexbins with 25km2 area -- these could be smaller but that would increase runtime 
cellsize <- 2 * sqrt(cellarea/((3*sqrt(3)/2))) * sqrt(3)/2 # calculate length across hexagon based on hexagon area
conus.grid <- st_make_grid(gridstart, cellsize = cellsize, square = F) %>% st_as_sf() %>% rename(geometry = x) %>% mutate(conus.grid.id = 1:nrow(.))

# create region
region <- make_region(rangelist = list(lake = GSL),
                      buffer = buffer,
                      grid = conus.grid)

ggplot() +
  geom_sf(data = region$sp.grid) +
  geom_sf(data = GSL, fill = "blue") +
  geom_sf(data = SLC, fill = "brown")


# Covariates ----

# download covariates
setwd("../species-futures/")
land <- get_landcover(locs = region$sp.grid, path = "data/USA/", id.label = "conus.grid.id")
footprint <- get_footprint(locs = region$sp.grid, id.label = "conus.grid.id")
elev <- get_elevation(locs = region$sp.grid, path = "data/USA/", id.label = "conus.grid.id")
setwd("../flexiSDM-great-salt-lake-cyclops/")

# join covariates together
covar <- full_join(land, footprint, by = "conus.grid.id") %>%
  full_join(elev, by = "conus.grid.id") %>%
  select(conus.grid.id, wetland, water, elevation, footprint, S) %>%
  rename(grid.id = conus.grid.id)

# There are a couple of hexbins with really high wetland values. To make this
# example more straightforward, truncate these values. This is not necessarily
# how you would handle this situation in a real analysis.
covar$wetland[which(covar$wetland > 0.5)] <- 0.5


covar_unscaled <- covar

# scale numeric cols
numcols <- sapply(covar, is.numeric)
numcols <- which(numcols)
covar[,numcols] <- sapply(covar[,numcols], scale_this)

summary(covar)

# visualize covariates
covlabs <- data.frame(covariate = c("wetland", "elevation", "footprint", "S"),
                      Label = c("Wetland", "Elevation", "Footprint", "South-facing"))

plot_covar(rename(covar, conus.grid.id = grid.id),
                  region,
                  cov.labs = covlabs,
                  scaled = T)$plot

cor_covar(rename(covar, conus.grid.id = grid.id), 
                 covlabs)$plot


# Simulate species ----

# Now we are going to create some random parameter values. Because this is
# a fake species, these values don't really matter. The point is that we
# know what the true values are, so we can compare our model estimates
# to the truth.

# Specify parameters
covs.z <- c("wetland", "footprint", "elevation", "elevation2", "S")
params <- c(0.46, -0.34, 0.27, -0.06, -0.21)
names(params) <- covs.z

# create quadratic term
covar <- covar %>%
  mutate(elevation2 = elevation * elevation) %>%
  select(all_of(c("grid.id", covs.z)))

# rename grid.id
grid <- region$sp.grid %>% rename(grid.id = conus.grid.id)

# create spatial noise
# we want this to be small compared to XB, so divide it by 10000
spat <- data.frame(grid.id = grid$grid.id,
                   spat = sim_cov(locs = grid, range = 100000)/10000)
summary(spat$spat)

# put together use list -- this tells sim_species() what values
# to use when simulating the species. Anything not included
# in `use` will be simulated within the function.
use <- list(grid.sim = grid, 
            grid.inf = grid,
            covar = covar,
            param = params,
            spat = spat)

# simulate species!
sp <- sim_species(use = use)


map_species_data(
  title = "Simulated species",
  region = region,
  plot = "true",
  sim = sp
)$plot

summary(sp$true$intensity)


# Simulate data ----

## iNat (PO) ----

# Create effort metric - slightly correlated with human footprint but with noise
effort <- (sim_cov(locs = grid, range = 10) + footprint$footprint/10)
ggplot(grid) + geom_sf(aes(fill = effort))

# create use list -- this tells sim_POdat() what values
# to use when simulating the PO records. Anything not included
# in `use` will be simulated within the function.
covar.eff <- data.frame(grid.id = grid$grid.id,
                        effort = scale(effort))
param <- c("effort" = 3.06)
use <- list(covar = covar.eff, param = param, alpha = 0.01)

# Simulate PO. This may take a couple minutes.
start <- Sys.time()
inat <- sim_POdat(
  true      = sp,
  use       = use,
  cell.area = 0.005          
)
Sys.time() - start
nrow(inat$dat) # number of PO records


## Northern Salt Lake Beast Brigade (DND) ----

# start with GSL, add buffer, then crop to above 41.2 lat
buffGSL <- GSL %>%
  st_transform(3857) %>% st_buffer(80000) %>% 
  st_transform(4326)
cropbb <- buffGSL %>% st_bbox()
cropbb["ymin"] <- 41.2
northGSL <- st_crop(buffGSL, cropbb)

# find grid cells in that region
northGSLgrid <- st_intersection(st_transform(region$sp.grid, 4326), northGSL) %>%
  select(conus.grid.id, geometry) %>%
  rename(grid.id = conus.grid.id)

# create locations within that region
locs <- sim_survey_sites(grid = st_transform(northGSLgrid, 3857), 
                         method = "random", nsite = 100, nvisit = 1)

# visualize
ggplot() +
  geom_sf(data = region$sp.grid, aes(fill = sp$true$intensity, color = sp$true$intensity)) +
  geom_sf(data = northGSLgrid, fill = NA, color = "white") +
  geom_sf(data = northGSL, fill = NA, color = "red") +
  geom_sf(data = locs, color = "green")

# create use list -- this tells sim_surveydat() what values
# to use when simulating the survey. Anything not included
# in `use` will be simulated within the function.
use <- list(locs = locs)

# simulate
NSLBB <- sim_surveydat(
  true = sp,
  use = use, 
  survey.type = "DND", # detection/nondetection survey
  
  # since detection covariates were not in the use list, the function
  # will simulate them internally using these instructions:
  covFormSite = "linear", # 1 linear site-level detection covariate
  covFormVisit = ""       # 0 visit-level detection covariates
)





## Utah Mythical Creature Count (count) ----

grid1 <- region$region %>% st_transform(crs = 4326)
cropbb <- grid1 %>% st_bbox()
cropbb["ymax"] <- 42 # survey ends at northern Utah border
grid1 <- st_crop(grid1, cropbb) %>%
  st_intersection(st_transform(region$sp.grid, 4326), region$sp.grid) %>%
  select(conus.grid.id, geometry) %>%
  rename(grid.id = conus.grid.id)

# simulate survey sites regularly spaced across Utah
locs <- sim_survey_sites(grid = st_transform(grid1, 3857), 
                         method = "regular", nsite = 50, nvisit = 1)

use <- list(locs = locs)

UMCC <- sim_surveydat(
  use = use,
  true = sp,
  survey.type = "count", # count survey
  
  # same detection covariates as DND
  covFormSite = "linear",
  covFormVisit = ""
)



## SLC Biodiversity Survey (count, multiple visits) ----

# start with SLC, add buffer
SLCgrid <-  SLC %>%
  st_transform(3857) %>% st_buffer(10000) %>% 
  st_transform(4326) %>% 
  st_intersection(st_transform(region$sp.grid, 4326)) %>%
  select(conus.grid.id, geometry) %>%
  rename(grid.id = conus.grid.id)

# simulate random survey sites within SLC
locs <- sim_survey_sites(grid = st_transform(SLCgrid, 3857), 
                         method = "random", nsite = 20, nvisit = 10)

use <- list(locs = locs)

# simulate
SLCBS <- sim_surveydat(
  true = sp,
  use = use,
  survey.type = "count", 
  covFormSite = "linear",
  covFormVisit = "", 
  link = "logit" # needs logit link because it uses an n-mixture style model
)



# Format data and visualize ----
datalist <- list(PO_iNat = inat,
                 DND_NSLBB = NSLBB,
                 count_UMCC = UMCC,
                 count_SLCBS = SLCBS)
sim_out <- format_sim(true = sp,
                      data = datalist)

region       <- sim_out$region
species.data <- sim_out$species.data
covar        <- sim_out$covar


pl <- map_species_data(title = "Simulated data",
                 region = region,
                 plot = "samples",
                 sim = T,
                 details = T,
                 species.data = species.data,
                 inat.agg = FALSE,
                 plot.range = FALSE,
                 plot.region = TRUE,
                 subtitle = FALSE)$plot
ggsave(pl, file = "data/datamap.jpg", height = 8, width = 8)

pl <- map_species_data(title = "True distribution",
                 region = region,
                 plot = "true",
                 sim = sp)$plot
ggsave(pl, file = "data/trueabundance.jpg", height = 8, width = 8)

true <- sp

save(region, datalist, species.data, covar, true, file = "data/data.RData")

