


# confirm correct package version - remove this before ESA
if (packageVersion("SpFut.flexiSDM") != "1.1.1") {
  remotes::install_git(
    "https://code.usgs.gov/eastern-ecological-science-center/nearmi/SpFut-flexiSDM.git",
    ref = "1.1.1"
  )
}


# This script fits a model for the Great Salt Lake Cyclops and compares
# the model output with the true distribution and parameters.

# ---------------------------------------------------------------------------
#
# SpFut.flexiSDM fits an integrated species distribution model (SDM): a single
# Bayesian model that combines several datasets of DIFFERENT types (e.g.,
# presence-only, detection/non-detection, counts) to estimate one shared
# distribution for a species. It can flexibly incorporation datasets with
# distinct observation models while all datasets inform the same
# underlying ecological process. The model is written in NIMBLE and fit with MCMC.
#
# Two layers to keep straight throughout this script:
#   - PROCESS model: the true ecology. 
#   - OBSERVATION model: how each dataset observes the process. 
#
# Because this is a SIMULATED species, we know the true parameters and the true
# distribution, so at the end we can check how well the model recovered them.
#
# Overall workflow:
#   1. Load and explore the core data inputs
#   2. Build model components
#   3. Assemble NIMBLE inputs
#   4. Fit the model and summarize the samples
#   5. Inspect chains and parameter estimates
#   6. Map estimated intensity and occupancy
#   7. Compare estimates to the known truth
#   8. Compare the integrated model to models fit to each dataset alone
# ---------------------------------------------------------------------------



# Set up environment ----
library(tidyverse)
library(SpFut.flexiSDM)

# `data.RData` holds the four components this workflow revolves around:
#   region        - the study area and its spatial grid (hexbins)
#   species.data  - the observation datasets (iNat, NSLBB, UMCC, SLCBS)
#   covar         - covariate values for every hexbin
#   true/datalist - the true simulated distribution and parameters (for validation)
load("data/data.RData", verbose = T)

# Fix the random seed so the model fit is reproducible
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

# There are two types of covariates used here:
#   - covs.z : PROCESS covariates. These drive the species' true distribution
#              and their coefficients are shared across all datasets.
#   - covs.PO: OBSERVATION covariates for the presence-only dataset,
#              e.g. sampling effort. These explain where people LOOKED, not
#              where the species IS.

# Define linear and quadratic covariates. These need to match the column names
# in `covar`. To include a quadratic term, add the squared column ("elevation2")
# alongside the linear one ("elevation").
covs.z    <- c("wetland", "footprint", "S", "elevation", "elevation2")

# Specify covariates for iNat data. Because PO data are modeled at the hexbin level,
# these covariates should also be stored in `covar`
covs.PO <- c("effort")


# covar is a data.frame that contains covariate information
#     - conus.grid.id must match conus.grid.id in region$species.grid

head(covar)

# NOTE: `covar` needs to contain scaled values and squared values for any
# quadratic terms


# Map some covariates
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

# Before NIMBLE can run, the spatial units need bookkeeping: an index that ties
# each hexbin to a row number, a train/test split for cross-validation, and
# (optionally) a coarser grid for the spatial random effect. That is what the
# next two sub-sections build.

## gridkey ----

# The gridkey contains:
#    - conus.grid.id, which maps each hexbin to the region$sp.grid sf object
#    - grid.id, which is used to index hexbins in the nimble code
#    - group, which defines each hexbin as "test" or "train"


# First, generate blocks to use for cross validation. Spatial (block) CV holds
# out whole regions of space rather than random cells, which gives a more
# honest test of how well the model predicts to unsampled areas. `rows`/`cols`
# lay a grid of blocks over the region; `k` is the number of folds those blocks
# are grouped into.
spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)

head(spatblocks)

ggplot() +  
  geom_sf(data = spatblocks, aes(fill = as.factor(folds))) + # CV blocks
  geom_sf(data = region$region, fill = NA) # overlay study region


# Assign a fold to leave out (options are "none", or 1:k).
# "none" fits the model to ALL cells (no hold-out) -- the usual choice for a
# final fit. Setting fold.out to 1, 2, ... instead reserves that fold as test
# data so you can evaluate out-of-sample prediction.
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

# Five things are fed into the NIMBLE model, built below:
#   data      - the observed values (counts, detections, etc.)
#   constants - fixed quantities the model needs (indices, covariate matrices,
#               dataset sizes) that are NOT being estimated
#   code      - the model definition itself (priors + likelihood), written
#               automatically by nimble_code() from `data` + `constants`
#   inits     - starting values for the MCMC chains
#   params    - which parameters to monitor (save) from the MCMC run



## Data ----

# Define `file.info`, which describes each dataset. This is the table that tells
# the model HOW to treat each dataset -- its observation type and which of its
# columns are detection covariates. One row per dataset in species.data$obs.
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


# Getting to NIMBLE-ready inputs is a two-step process:
#   1. sppdata_for_nimble() reshapes EACH observation dataset to prep it for NIMBLE.
#   2. data_for_nimble() then bolts those onto the shared PROCESS model and
#      produces the final `data` + `constants`.

# Step 1: format each dataset in species.data for nimble.
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

# sp.data now contains `data` and `constants` for each dataset. The names encode
# the data type and position: PO1 = presence-only (iNat), DND2 = detection/
# non-detection (NSLBB), count3/count4 = the two count datasets (UMCC, SLCBS).
sp.data$PO1
sp.data$DND2
sp.data$count3
sp.data$count4


# Step 2: combine the per-dataset observation pieces with the shared PROCESS
# (distribution) model
tmp <- data_for_nimble(sp.data = sp.data,
                       region = region,
                       gridkey = gridkey,
                       covar = covar,
                       covs.z = covs.z,
                       
                       # spatial model
                       sp.auto = TRUE,      
                       coarse.grid = FALSE,
                       spatRegion = NULL)

data      <- tmp$data
constants <- tmp$constants

# `data` and `constants` are ready for nimble
data
constants


## NIMBLE parameters ----

# `nimble_code()` writes the NIMBLE model code for you,
# tailored to the datasets and structure implied by `data` and `constants`.
code <- nimble_code(data = data,
                    constants = constants,
                    
                    # other inputs
                    Bprior = "dnorm(0,1)",
                    block.out = fold.out,
                    
                    # what types of observation models to use
                    min.visits.incl = 3,
                    occ.mod = TRUE,
                    nmix.mod = TRUE,
                    
                    # spatial model
                    sp.auto = TRUE,
                    coarse.grid = FALSE,
                    zero_mean = TRUE,
                    tau = 1,
                    
                    # where to write nimble code to
                    path = tempdir())

# `data` and `constants` are also used to generate the initial values and
# a vector of the parameters to save.

# inits is a FUNCTION of a seed (x): each MCMC chain calls it with a different
# seed to get its own randomized starting values.
inits <- function(x) {nimble_inits(data = data,
                                   constants = constants,
                                   sp.auto = TRUE,
                                   min.visits.incl = 3,
                                   occ.mod = TRUE,
                                   nmix.mod = TRUE,
                                   seed = x)}

# params lists which quantities to monitor. lambda = TRUE and XB = TRUE also
# save the derived per-cell intensity and linear predictor, which we map later.
params <- nimble_params(data = data,
                        constants = constants,
                        lambda = TRUE,
                        XB = TRUE,
                        sp.auto = TRUE,
                        effort = FALSE)


# Fit the model ----

## Fit the model ----

# MCMC settings:
#   iter   - total iterations per chain
#   thin   - keep every xth sample (reduces autocorrelation/storage)
#   burnin - discard the first 75% as warm-up before the chain has converged
# These are small for a quick demo; a real analysis typically needs many more.
iter   <- 5000
thin   <- 5
burnin <- floor(iter * 0.75)

# Fitting takes a while. Set to FALSE to skip it and load saved output instead
# (see the `if (!exists("out"))` block below).
run.during.break <- TRUE

if (run.during.break) {

  # nimbleParallel() runs the MCMC chains in parallel and returns the posterior
  # samples. This is the actual model-fitting step.
  samples <- nimbleParallel(code = code,
                            data = data,
                            constants = constants,
                            inits = inits,
                            param = params,
                            iter = iter,
                            burnin = burnin,
                            thin = thin)
  
  
  ## Summarize model output ----

  # get_derived() adds derived quantities (e.g. per-cell occupancy)
  # to the raw samples of each chain, computed from the monitored parameters.
  samples <- lapply(samples,
                    get_derived,
                    data = data,
                    
                    # projections (see vignette)
                    project = 0,
                    proj.data = data,
                    
                    # spatial effect
                    sp.auto = TRUE,
                    coarse.grid = FALSE,
                    spatRegion = spatRegion)

  # summarize_samples() collapses the chains into tidy per-parameter summaries
  # (posterior means, credible intervals, Rhat, ESS). `out` is what the plotting
  # and comparison functions below consume.
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


# If you were unable to fit the model (e.g. run.during.break = FALSE, or the fit
# failed), use a pre-computed `out` saved on disk so the rest of the
# script still runs.
if (!exists("out")) {
  load("data/output.RData")
}


# View model output ----

# Two complementary diagnostics are used throughout this section:
#   - plot_chains(): the raw MCMC traceplots. Use these to check convergence --
#     the chains should look like fuzzy caterpillars overlapping each other.
#   - plot_pars():   the summarized posterior (point estimate + credible
#     interval) for each parameter.
# Recall the parameter groups: "B" = process coefficients, "alpha" = dataset
# intercepts, "observation" = detection parameters.

# view beta (process) parameters -- the shared covariate effects on distribution
plot_chains(samples, data = data, cov.labs = cov.labs,
                       plot = "B", cutoff = 0)$plot
plot_pars(out = out$process.coef, cov.labs = cov.labs)$plot
# plot_effects() shows the fitted response curve for each covariate, on the
# original (unscaled) covariate axis.
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

# map_species_data() draws the model's spatial predictions across the hexbins.
# `plot = "lambda"` maps relative intensity/abundance; `plot = "psi"` maps
# occupancy probability. sim = TRUE signals that this is a simulated species.

# Plot relative abundance
map_species_data(title = "Estimated intensity",
                 region = region,
                 out = out,
                 sim = TRUE,
                 plot = "lambda",
                 plot.range  = FALSE,
                 plot.region = TRUE)$plot

# Plot relative occupancy probability
map_species_data(title = "Estimated occupancy probability",
                 region = region,
                 out = out,
                 sim = TRUE,
                 plot = "psi",
                 plot.range = FALSE,
                 plot.region = TRUE)$plot

# See documentation for many more options to map
?map_species_data()

# Compare model estimates to truth ----

# Because this species is simulated, we know the real answer and can grade the
# model. `sim_compare()` lines up the estimates in `out` against the truth.
# It takes different "true" objects depending on what you're comparing:
#   - process / lambda / spat : compare against `true` (the process truth)
#   - alpha / obs             : compare against `datalist` (per-dataset truth)
# The `plot` argument picks which quantity to compare.

# datalist stores the true observation-model values per dataset, but in the same
# order as species.data$obs -- so copy those names over before using it.
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

# Does combining all datasets beat fitting each dataset on its own? 
# The saved objects below hold sim_compare() results from single-dataset models 
# (one per dataset). We append the results from our integrated ("all") model 
# and plot them side by side against truth.

# load outputs for models that fit each dataset individually.
# Provides: all.proc, all.lamb, all.obs, all.alph, all.spat
load("data/model-outputs.RData", verbose = T)

# each of these objects contains sim_compare() output for each model


## Compare process parameters ----
# Pull the integrated model's process estimates and tag them model = "all"
proc <- sim_compare(out, true = true, plot = "process")$dat %>%
  mutate(model = "all") # label this model

# then stack them onto the single-dataset results for a combined comparison.
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

# Same idea as above, but comparing the per-cell intensity (lambda) predictions.
# In the scatterplot, each point is a hexbin (estimated vs. true intensity);
# points falling on the 1-to-1 line mean the model correctly estimated the
# value in that cell.
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
q90 <- quantile(all.lamb.sf$mean, 0.90)
all.lamb.sf$mean[which(all.lamb.sf$mean > q90)] <- q90

# map
ggplot(all.lamb.sf) +
  geom_sf(aes(fill = mean, color = mean)) +
  facet_wrap(~ model) +
  
  # color scale
  viridis::scale_fill_viridis(option = "magma") +
  viridis::scale_color_viridis(option = "magma")


# map differences (estimate - truth) to see WHERE each model is off.
ggplot(all.lamb.sf) +
  geom_sf(aes(fill = mean - true, color = mean - true)) +
  facet_wrap(~ model) +

  # color scale
  scale_fill_gradient2() +
  scale_color_gradient2()
