


# confirm correct package version - remove this before ESA
if (packageVersion("SpFut.flexiSDM") != "1.1.1") {
  remotes::install_git(
    "https://code.usgs.gov/eastern-ecological-science-center/nearmi/SpFut-flexiSDM.git",
    ref = "1.1.1"
  )
}


# This script fits a model for the Great Salt Lake Cyclops and compares
# the model output with the true distribution and parameters.



# Set up environment ----
library(tidyverse)
library(SpFut.flexiSDM)

load("data/data.RData", verbose = T)

set.seed(1)


# Explore inputs ----

## Region ----

# region is created by:
#     - SpFut.flexiSDM::make_region() for real species
#     - SpFut.flexiSDM::format_sim() for simulated species


# It has four components:

# For a simulated species, range is NULL
# For a real species, range is the boundary that was used to create the study region
region$range


# region is the outer boundary of the study region
ggplot(region$region) + geom_sf()


# sp.grid contains the spatial units for analysis
# This example uses hexbins with area 25km2
ggplot(region$sp.grid) + geom_sf()
nrow(region$sp.grid)
# conus.grid.id is a unique label for each hexbin
head(region$sp.grid)



# boundary dictates where the region is cut off
# In this case, it is the same as region, but it can be different
# E.g., you could use a state boundary to restrict the study region to within the state
ggplot(region$boundary) + geom_sf()


## Species data ----

# species.data is created by:
#     - SpFut.flexiSDM::load_species_data() for real species
#     - SpFut.flexiSDM::format_sim() for simulated species


# It has two elements:

# locs: the locations of each data point in discrete (disc) and continuous (cont) space
#       - locs$disc (data.frame) contains each point assigned to a hexbin
#       - locs$cont (sf) contains each point assigned to coordinates

head(species.data$locs$disc)

head(species.data$locs$cont)

ggplot() +
  geom_sf(data = region$region) +  # start with region
  geom_sf(data = species.data$locs$cont, aes(color = source)) # overlay points


# obs: list of elements that each contain data from one dataset
#      - each element contains date, location, and species count info for one dataset

head(species.data$obs$iNat)
head(species.data$obs$NSLBB)
head(species.data$obs$UMCC)
head(species.data$obs$SLCBS)


# Use SpFut.flexiSDM::map_species_data() with region and species.data
map_species_data(title = "Simulated data",
                 region = region,
                 plot = "samples",
                 sim = TRUE,
                 details = TRUE,
                 species.data = species.data,
                 inat.agg = FALSE,
                 plot.range = FALSE,
                 plot.region = TRUE,
                 subtitle = FALSE)$plot


## Covariates ----

# Define linear and quadratic covariates. These need to match the column names
# in `covar`
covs.z    <- c("wetland", "footprint", "S", "elevation", "elevation2")

# Specify covariates for iNat data. Because PO data are modeled at the hexbin level,
# these covariates should also be stored in `covar`
covs.PO <- c("effort")


# covar is a data.frame that contains covariate information
#     - conus.grid.id must match conus.grid.id in region$species.grid

head(covar)

# NOTE: `covar` needs to contain scaled values and squared values for any
# quadratic terms



gridcovar <- full_join(region$sp.grid, covar, by = "conus.grid.id")

ggplot(gridcovar) +
  geom_sf(aes(fill = wetland))

ggplot(gridcovar) +
  geom_sf(aes(fill = elevation))



# cov.labs assigns labels to each covariate
cov.labs <- data.frame(covariate = names(covar),
                       Label = c("conus.grid.id", "Wetland", "Footprint", "Elevation",
                                 "Elevation^2", "South", "Effort"))  %>%
  filter(!covariate == "conus.grid.id")
cov.labs


# Plot covariates
plot_covar(covar = covar, region = region, cov.labs = cov.labs, scaled = TRUE)$plot
cor_covar(covar = covar, cov.labs = cov.labs, color.threshold = 0.25)$plot





# Set up model components ----

## gridkey ----

# The gridkey contains:
#    - conus.grid.id, which maps each hexbin to the region$sp.grid sf object
#    - grid.id, which is used to index hexbins in the nimble code
#    - group, which defines each hexbin as "test" or "train"


# First, generate blocks to use for cross validation. `k` is the number
# of folds that blocks are grouped into
spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)

head(spatblocks)

ggplot() +  
  geom_sf(data = spatblocks, aes(fill = as.factor(folds))) + # CV blocks
  geom_sf(data = region$region, fill = NA) # overlay study region


# Assign a fold to leave out (options: none, 1:k)
fold.out <- "none"

# Now generate gridkey, assigning which fold to leave out
gridkey <- make_gridkey(region, spatblocks, fold.out = fold.out)

head(gridkey)

key <- full_join(region$sp.grid, gridkey, by = "conus.grid.id")

ggplot() +
  geom_sf(data = key, aes(fill = group))


# If you use `fold.out = 1`, the cells in fold 1 will become test data
gridkey1 <- make_gridkey(region, spatblocks, fold.out = 1)

head(gridkey1)

key1 <- full_join(region$sp.grid, gridkey1, by = "conus.grid.id")

ggplot() +
  geom_sf(data = key1, aes(fill = group))


## spatRegion ----

# To improve processing time, we can create a coarse spatial grid that
# we use to fit the ICAR model

spatRegion <- suppressWarnings(make_spatkey(region$sp.grid))

# spatRegion contains two elements:


# spatkey: assigns each conus.grid.id to a spat.grid.id
head(spatRegion$spatkey)

# spatgrid: sf object containing spatial grid
head(spatRegion$spat.grid)

ggplot(spatRegion$spat.grid) + geom_sf()
# rather than estimating a spatial random effect for each hexbin (3213),
# using the coarse grid, random effects are only estimated for each group of hexbins (490)

# Alternatively, set spatRegion = NULL to estimate a random effect for each hexbin
spatRegion <- NULL



# Set up NIMBLE model ----

## Data ----

# Define `file.info`, which describes each dataset.
# `file.info` must contain:

#   - `file.label`: the label that should be used for the dataset; must match
#     names(species.data$obs)

#   - `covar.mean`: comma-separated detection covariate(s) that should be
#     averaged across passes; must match column name in species.data$obs$<file.label> file

#   - `covar.sum`: comma-separated detection covariate(s) that should be
#     summed across passes; must match column name in data file

#   - `data.type`: data type of dataset; must be "PO", "DND", "or "count"

#   - `PO.extent`: describes the spatial extent of PO datasets; must be
#     "CONUS" or a two-letter state abbreviation; NA for non-PO datasets

file.info <- data.frame(file.label = c("iNat", "NSLBB", "UMCC", "SLCBS"),
                        data.type = c("PO", "DND", "count", "count"),
                        covar.mean = c("effort", "cov.siteA", "cov.siteA", "cov.siteA"),
                        covar.sum = c(NA, NA, NA, NA),
                        PO.extent = c("CONUS", NA, NA, NA))
file.info


# format species.data for nimble
sp.data <- sppdata_for_nimble(species.data = species.data,
                              region = region,
                              covar = covar,
                              
                              # some things that were defined earlier
                              file.info = file.info, 
                              covs.PO = covs.PO,
                              keep.conus.grid.id = gridkey$conus.grid.id[gridkey$group == "train"],

                              # what types of observation models to use
                              occ.mod = TRUE,
                              nmix.mod = TRUE,
                              min.visits.incl = 3,
                              
                              # these options only apply to real data; see vignette
                              stategrid = NULL,    # no state grid for simulated data
                              statelines.rm = FALSE)

# sp.data now contains `data` and `constants` for each dataset
sp.data$PO1
sp.data$DND2
sp.data$count3
sp.data$count4


# Now combine that species data with distribution process data
tmp <- data_for_nimble(sp.data = sp.data,
                       covar = covar,
                       covs.z = covs.z,
                       region = region,
                       gridkey = gridkey,
                       sp.auto = TRUE,      
                       coarse.grid = FALSE,
                       spatRegion = NULL)

data      <- tmp$data
constants <- tmp$constants

# `data` and `constants` are ready for nimble
data
constants


## NIMBLE parameters ----

# `nimble_code()` writes the NIMBLE code for the model based
# on the information in `data` and `constants`.

code <- nimble_code(data = data,
                    constants = constants,
                    path = tempdir(),
                    sp.auto = TRUE,
                    coarse.grid = FALSE,
                    Bprior = "dnorm(0,1)",
                    block.out = fold.out,
                    zero_mean = TRUE,
                    tau = 1,
                    min.visits.incl = 3,
                    occ.mod = TRUE,
                    nmix.mod = TRUE)

# `data` and `constants` are also used to generate the initial values and
# a vectors of the parameters to save.

inits <- function(x) {nimble_inits(data = data,
                                   constants = constants,
                                   sp.auto = TRUE,
                                   min.visits.incl = 3,
                                   occ.mod = TRUE,
                                   nmix.mod = TRUE,
                                   seed = x)}

params <- nimble_params(data = data,
                        constants = constants,
                        lambda = TRUE,
                        XB = TRUE,
                        sp.auto = TRUE,
                        effort = FALSE)


# Fit the model ----

## Fit the model ----

# define model parameters
iter   <- 5000
thin   <- 5
burnin <- floor(iter * 0.75)

run.during.break <- TRUE

if (run.during.break) {
  
  samples <- nimbleParallel(code = code,
                            data = data,
                            constants = constants,
                            inits = inits,
                            param = params,
                            iter = iter,
                            burnin = burnin,
                            thin = thin)
  
  
  ## Summarize model output ----
  samples <- lapply(samples, 
                    get_derived, 
                    data = data, 
                    project = 0,
                    proj.data = data,
                    sp.auto = TRUE,
                    coarse.grid = FALSE, 
                    spatRegion = spatRegion)
  
  
  out <- summarize_samples(samples = samples,
                           data = data,
                           constants = constants,
                           project = 0,
                           coarse.grid = FALSE,
                           block.out = fold.out,
                           gridkey = gridkey,
                           effort = FALSE,
                           cores = 2L)
  
}


# If you were unable to fit the model, read in the output:
if (!exists("out")) {
  load("data/output.RData")
}


# View model output ----

# view beta parameters
plot_chains(samples, data = data, cov.labs = cov.labs,
                       plot = "B", cutoff = 0)$plot
plot_pars(out = out$process.coef, cov.labs = cov.labs)$plot
plot_effects(data = data, out = out, cov.labs = cov.labs, unscale_covar = TRUE)$plot

# View chains for dataset intercepts (alpha)
plot_chains(samples, data = data, constants = constants, cov.labs = cov.labs,
                       plot = "alpha", cutoff = 0)$plot
plot_pars(out = out$alpha, cov.labs = cov.labs)$plot

# View chains for detection parameters
plot_chains(samples, data = data, constants = constants, cov.labs = cov.labs,
            plot = "observation", cutoff = 0)$plot
plot_pars(out = out$obs.coef, cov.labs = cov.labs)$plot


# Map estimates ----

# Plot relative abundance
map_species_data(title = "Estimated intensity",
                 region = region,
                 out = out,
                 sim = T,
                 plot = "lambda",
                 plot.range  = FALSE,
                 plot.region = TRUE)$plot

# Plot relative occupancy probability
map_species_data(title = "Estimated occupancy probability",
                 region = region,
                 out = out,
                 sim = T,
                 plot = "psi",
                 plot.range = FALSE,
                 plot.region = TRUE)$plot

# See documentation for many more options to map
?map_species_data()

# Compare model estimates to truth ----

# Use `sim_compare()` to compare the true parameters with the parameters
# estimated by the model.


names(datalist) <- names(species.data$obs)

# sim_compare() (and most flexiSDM plotting functions) produce a list with 
# two objects:
#     - plot: a ggplot object
#     - data: a data.frame containing the data used to make the plot
tmp <- sim_compare(out, true = true, plot = "process")

tmp$plot
tmp$data

# look at outputs
sim_compare(out, true = true, plot = "process")$plot
sim_compare(out, true = true, plot = "lambda")$plot
sim_compare(out, true = true, plot = "spat")$plot

sim_compare(out, true = datalist, plot = "alpha")$plot
sim_compare(out, true = datalist, plot = "obs")$plot



# Compare model estimates with SDMs fit with each individual dataset ----

# load outputs for models that fit each dataset individually
load("data/model-outputs.RData", verbose = T)

# each of these objects contains sim_compare() output for each model


## Compare process parameters ----
proc <- sim_compare(out, true = true, plot = "process")$dat %>%
  mutate(model = "all") # label this model

# add to all.proc
all.proc <- bind_rows(all.proc, proc)


ggplot(all.proc) +
  # plot 0
  geom_hline(yintercept = 0) +
  
  # plot true value
  geom_hline(aes(yintercept = true), color = "darkblue", linetype = "dashed", linewidth = 0.75) +
  
  # plot estimates with CI
  geom_pointrange(aes(x = covariate, y = mean, ymin = lo, ymax = hi, 
                      color = model, group = model), 
                  position = position_dodge(width = 0.8)) +
  
  # format
  facet_wrap(~ covariate, scales = "free_x") +
  theme_bw()




## Compare intensity estimates ----

lamb <- sim_compare(out, true = true, plot = "lambda")$dat %>%
  mutate(model = "all")
all.lamb <- bind_rows(all.lamb, lamb)

ggplot(all.lamb) +
  
  # plot lambdas
  geom_point(aes(x = true, y = mean, color = model)) +
  
  # 1-to-1 line
  geom_abline(intercept = 0, slope = 1) +
  
  # parameters
  facet_wrap(~ model, scales = "free") +
  coord_cartesian(xlim = c(0, 25), ylim = c(0, 25))



# map estimated intensity 

# join with sp.grid to make sf object
all.lamb.sf <- region$sp.grid %>%
  full_join(all.lamb, by = "conus.grid.id")

# truncate at 95%ile to help color scale
q95 <- quantile(all.lamb.sf$mean, 0.90)
all.lamb.sf$mean[which(all.lamb.sf$mean > q95)] <- q95

# map
ggplot(all.lamb.sf) +
  geom_sf(aes(fill = mean, color = mean)) +
  facet_wrap(~ model) +
  
  # color scale
  viridis::scale_fill_viridis(option = "magma") +
  viridis::scale_color_viridis(option = "magma")


# map differences
ggplot(all.lamb.sf) +
  geom_sf(aes(fill = mean - true, color = mean - true)) +
  facet_wrap(~ model) +
  
  # color scale
  scale_fill_gradient2() +
  scale_color_gradient2()
