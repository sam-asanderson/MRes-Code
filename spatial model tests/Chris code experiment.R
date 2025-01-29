rm(list=ls())

## Set seed
set.seed(1)

####################################
## Parameters
####################################
## Simulation length (days)
t_sim <- 100

## Size of the landscape (x*x)
size <- 10

## Number of pats per day
pats_per_day <- 2

## Proportion of pats treated
n_treat <- 0.5

## Number of species (placeholder)
n_spp <- 2

## Starting abundances
start_abun <- c(spp_1 = 500, spp_2 = 500)

## Generation times (days)
spp_gen <- c(spp_1 = 60, spp_2 = 180)

## Resources of a pat
resources <- 500

## Lifetimes of pats
treated_lifetime <- 180
untreated_lifetime <- 90

####################################
## Object setup
####################################

## Set up empty field
field <- matrix(0, nrow = size, ncol = size)

## Empty object to store locations of pats
pats_loc <- array(0, dim = c(size, size, t_sim))

## Empty object to store resources left in each pat
pats_res <- array(0, dim = c(size, size, t_sim))

## Treatment locations (ivermectin)
treat_loc <- array(0, dim = c(size, size, t_sim))

## Time left for each pat
pat_lifetime_left <- array(0, dim = c(size, size, t_sim))

## Species abundances at locations
## First dimension: x field
## Second: y field
## Third: time
## Fourth: species
spp_abun <- array(0, dim = c(size, size, t_sim, n_spp))

####################################
## Starting conditions
####################################

for (l in 1:n_spp) {
  ## Evenly distribute initial abundances across the grid
  spp_abun[,,1,l] <- start_abun[l] / length(spp_abun[,,1,l])
}

####################################
## Simulation
####################################

for (i in 1:t_sim) {
  ## Reduce all pat lifetimes by 1 each time step, but not below zero
  if (i > 1) { 
    pat_lifetime_left[,,i] <- pmax(pat_lifetime_left[,,i - 1] - 1, 0)
  }
  
  ##############################################
  ## Step 1 - Add new pats
  pats_current <- pats_loc[,,i]
  
  ## Randomly select locations for new pats
  pats_t_I <- sample(size * size, pats_per_day)
  
  ## Add new pats to the field
  pats_current[pats_t_I] <- 1
  
  ## Update pats_loc for the current time step
  pats_loc[,,i] <- pats_current
  
  ##############################################
  ## Step 2 - Assign treatment to new pats
  
  ## Pull out the current treatment data
  treat_current <- treat_loc[,,i]
  
  ## Determine which new pats are treated (1) or untreated (0)
  treat_pats <- rbinom(pats_per_day, 1, n_treat)
  
  ## Assign treatments to the new pats
  treat_current[pats_t_I] <- treat_pats
  
  ## Update treatment array
  treat_loc[,,i] <- treat_current
  
  ## Assign lifetimes to new pats
  pat_lifetime_left[,,i][pats_t_I] <- ifelse(treat_current[pats_t_I] == 1, treated_lifetime, untreated_lifetime)
  
  ##############################################
  ## Step 3 - Species simulations
  
  # spp_abun[,,i,l] 
  
}

## Output recorded data
print(pat_lifetime_left[,,1:10])

