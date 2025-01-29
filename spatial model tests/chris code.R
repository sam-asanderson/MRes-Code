rm(list=ls())

##set seed
set.seed(1)

####################################
##parameters
####################################
##sim length (days)
t_sim <- 100

##size of the landscape (x*x)
size <- 10

##number of pats per day
pats_per_day <- 2

##proportion of pats treated
n_treat <- 0.5

##n_spp - placeholder for the moment to create array
n_spp <-2

##starting abundances
start_abun <- c(spp_1=500, spp_2 = 500)

##generation times (days)
spp_gen <- c(spp_1=60, spp_2 = 180)
##resources of a pat
resources <- 500

##lifetimes of pats
treated_lifetime = 180
untreated_lifetime = 90
####################################
##object set up
####################################

##set up empty field
field <- matrix(0, nrow=size, ncol=size)

##empty object to store locations of pats
pats_loc <- array(0, dim=c(size, size, t_sim))

##empty object to store resources left in each pat
pats_res <- array(0, dim=c(size, size, t_sim))

##treatment locations (ivermectin)
treat_loc <- array(0, dim=c(size, size, t_sim))

##time left
pat_lifetime_left <- array(0, dim=c(size, size, t_sim))

##species abundances at locations
##first dimension is x field,?
##second is y field?
##third is time
##4th is species
spp_abun <- array(0, dim=c(size, size, t_sim, n_spp))

####################################
##starting conditions
####################################

for(l in 1:length(n_spp)){
  #l=1
  ##divide the starting abundance by the number of patches
  ##evenly distribute
  spp_abun[,,1,l]<-start_abun[l]/length(spp_abun[,,1,l])
}

####################################
##sim
####################################

##run sim
for(i in 1:length(t_sim)){
  ##define i
  i=1
  
  ##reduce all pat lifetimes by 1 each time step - ADD DON'T GO BELOW ZERO - EVENTUALLY WHEN THIS GETS TO ZERO, PATS WILL DISAPPEAR. 
  if (i > 1) { ##doesn't happen first time
    pat_lifetime_left[,,i] <- pat_lifetime_left[,,i - 1] - 1
  }
  ##############################################
  ##first step - add some pats
  pats_current <- pats_loc[,,i] 
  
  ##the locations of the pats being added this time
  pats_t_I <- sample(size*size, pats_per_day)
  
  ##add these pats to the field - 
  pats_current[pats_t_I] <-1
  
  pats_loc[,,i] <- pats_current
  ##############################################
  ##second step - define the infected pats
  
  ##pull out the current treatments
  treat_current <- treat_loc[,,i]
  
  ##define which of the new pats added today have ivermectin
  treat_pats <- rbinom(pats_per_day, 1, n_treat) #as two pats are added a day, two integers are recorded - either 1 or 0.  
  treat_current[pats_t_I] <- treat_pats #this updates treat_current with the 0 or 1 from treat_pats in the location given by pats_t_I of the new pats added that time step. 
  
  #pats at locations pats_t_I get updates in pat_lifetimes_left to 90 or 180 depending.
  pat_lifetime_left[pats_t_I] <- ifelse(treat_current[pats_t_I] == 1, treated_lifetime, untreated_lifetime) 
  
 
  ##############################################
  ##species simulations
  
  ##
  spp_abun[,,i,l]
  
  
}

print(pats_loc[,,1])
print(pat_lifetime_left[,,1:2])
dung_fly <- function()
print(pats_t_I)  
print(treat_current)
print(pats_current)
