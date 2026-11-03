
# This tutorial uses SpFut.flexiSDMv1.1.2. Using it with
# a different version may cause errors. For instructions on how
# to install SpFut.flexiSDMv1.1.2, see 00-set-up.R



# This script fits a model for the Great Salt Lake Cyclops using each
# individual dataset.


# Set up environment ----
library(tidyverse)
library(SpFut.flexiSDM)

load("data/data.RData")

set.seed(1)

# define covariates 
covs.PO <- c("effort") # covariate for PO iNat data
covs.z    <- c("wetland", "footprint", "S", "elevation", "elevation2") # process covariates

# Create spatial blocks
spatblocks <- make_CV_blocks(region, rows = 5, cols = 5, k = 3)

# Assign a fold to leave out (options: none, 1:k)
fold.out <- "none"

# Now, generate gridkey, assigning which fold to leave out 
gridkey <- make_gridkey(region, spatblocks, fold.out = fold.out)


# Create summary of each dataset - data type, names of covariates, and spatial extent
file.info <- data.frame(file.label = c("iNat", "NSLBB", "UMCC", "SLCBS"),
                        data.type = c("PO", "DND", "count", "count"),
                        covar.mean = c("effort", "cov.siteA", "cov.siteA", "cov.siteA"),
                        covar.sum = c(NA, NA, NA, NA),
                        PO.extent = c("CONUS", NA, NA, NA))


# Fit models ----
all.proc <- c()
all.alph <- c()
all.lamb <- c()
all.obs <- c()
all.spat <- c()
for (i in 1:4) {
  
  # pull data and set model name
  if (i == 1) {
    datalist1 <- list(PO_iNat = datalist$PO_iNat)
    model <- "iNat"
  } else if (i == 2) {
    datalist1 <- list(DND_NSLBB = datalist$DND_NSLBB)
    model <- "NSLBB"
  } else if (i == 3) {
    datalist1 <- list(count_UMCC = datalist$count_UMCC)
    model <- "UMCC"
  } else if (i == 4) {
    datalist1 <- list(count_SLCBS = datalist$count_SLCBS)
    model <- "SLCBS"
  }
  
  sim_out <- format_sim(true = true,
                        data = datalist1)
  
  region       <- sim_out$region
  species.data <- sim_out$species.data
  covar        <- sim_out$covar
  
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
  
  
  # Now combine that species data with distribution process data
  tmp <- data_for_nimble(sp.data = sp.data,
                         covar = covar,
                         covs.z = covs.z,
                         sp.auto = TRUE,      
                         coarse.grid = FALSE,
                         region = region,
                         gridkey = gridkey,
                         spatRegion = NULL)
  
  data      <- tmp$data
  constants <- tmp$constants
  
  
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
  
  
  ## Fit the model ----
  
  # define model parameters
  iter   <- 5000
  thin   <- 5
  burnin <- floor(iter * 0.75)
  
  
  
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
  
  
  names(datalist1) <- names(species.data$obs)
  
  # To bypass a bug in v1.1.2, remove covForm before
  # running sim_compare()
  true1 <- true
  true1$covForm <- NULL
  
  proc <- sim_compare(out, true = true1, plot = "process")$dat %>%
    mutate(model = model)
  
  lamb <- sim_compare(out, true = true1, plot = "lambda")$dat %>%
    mutate(model = model)
  
  spat <- sim_compare(out, true = true1, plot = "spat")$dat %>%
    mutate(model = model)
  
  alph <- sim_compare(out, true = datalist1, plot = "alpha")$dat %>%
    mutate(model = model)
  
  obs <- sim_compare(out, true = datalist1, plot = "obs")$dat %>%
    mutate(model = model)
  
  
  # save
  all.proc <- bind_rows(all.proc, proc)
  all.alph <- bind_rows(all.alph, alph)
  all.lamb <- bind_rows(all.lamb, lamb)
  all.obs <- bind_rows(all.obs, obs)
  all.spat <- bind_rows(all.spat, spat)
}


save(all.proc, all.alph, all.lamb, all.obs, all.spat, file = "data/model-outputs.RData")

