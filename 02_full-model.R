################################################################################
#' Simulation of Dung Insect Guild Population Dynamics Under Landscape
#' Management Scenarios
#'
#' Author: Sam Appleyard-Sanderson
#'
#' Date: 2025-07-21
################################################################################

# --- Load packages ------------------------------------------------------------
library(foreach)
library(doParallel)

# --- Clear environment & load functions ---------------------------------------
rm(list = ls()) # clear environment

source("2025-05-22 connectivity function.R")

# --- Setup parallel cluster ---------------------------------------------------
# detect how many cores are available (-1 for system)
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)
# using doParallelSNOW - name of backend used by foreach.
# using version 1.0.17.

# --- Set global parameters ----------------------------------------------------
t_sim <- 1000 # length of simulation (days)
size <- 44    # 100-600 dimension of pasture (no. patches = size * size)
guild_no <- 3 # no. of guilds
start_time <- 5 # no. of time steps that new populations are added
# creates connectivity values between all pats
connectivity_values <- GenerateConn(size = size, rate = 1, norm = TRUE)
# days spent in larval feeding stage
feeding_time <- c(
  # fly
  5,
  # dweller
  20,
  # tunneler
  20
)
# time spent pupating
pupation_time <- c(
  # fly
  15,
  # dweller
  20,
  # tunneler
  20
)
# eggs laid per female
fecundity <- matrix(c(
  # fly untreated
  35,
  # fly treated
  0,
  # dweller untreated
  15,
  # dweller treated
  0,
  # tunneler untreated
  15,
  # tunneler treated
  0
),
nrow = 3,
ncol = 2,
byrow = TRUE)
# biomass removed / individual / day
dung_removal <- matrix(
  c(
    # fly untreated
    0.03,
    # fly treated
    0.03,
    # dweller untreated
    0.1,
    # dweller treated
    0.05,
    # tunneler untreated
    0.25,
    # tunneler treated
    0.125
  ),
  nrow = 3,
  ncol = 2,
  byrow = TRUE
)

# --- Set treatment variables --------------------------------------------------
# no. new pats added / day
pats_per_day <- c(
  # SD 5, 10, 15 AU / hectare
  5, # 200
  10, # 400
  15 #600
)

tst_level <- seq(0, 1, 0.1) # prop. untreated pats
treatment_grid <-
  expand.grid(pats_per_day = pats_per_day, tst_level = tst_level)
n_repeats <- 1 # 10 # repeats per treatment
full_grid <- expand.grid(treatment_id = seq_len(nrow(treatment_grid)),
                         repeat_id = 1:n_repeats)

# --- Parallel simulation ------------------------------------------------------
# export custom functions and data to the cluster
clusterExport(cl, varlist = c("GenerateConn", "full_grid"))

# start of simulation
treatment_results <- foreach(tt = seq_len(nrow(full_grid)), .combine = rbind) %dopar% {
    t <- full_grid$treatment_id[tt]
    r <- full_grid$repeat_id[tt]
    pats_per_day <- treatment_grid$pats_per_day[t]
    tst_level <- treatment_grid$tst_level[t]
    cell <- list(
      pat_present = FALSE,
      pat_age = NA_integer_,
      pat_treatment = FALSE,
      adult_pat_abun = matrix(0, ncol = guild_no, nrow = 40),
      adult_field_abun = c(0, 0, 0),
      larvae_feeding_abun = matrix(0, ncol = guild_no, nrow = 20),
      larvae_pupating_abun = matrix(0, ncol = guild_no, nrow = 20),
      pat_biomass = 0
    )
    pasture <- matrix(
      # creates grid of cells
      data = lapply(1:(size * size), function(x)
        cell),
      nrow = size,
      ncol = size
    )
    pasture_results <- data.frame(
      # results dataframe updates per time step
      t_sim = numeric(t_sim),
      tst_level = numeric(t_sim),
      pats_per_day = numeric(t_sim),
      total_pat_number = numeric(t_sim),
      average_pat_biomass = numeric(t_sim),
      total_adult_pat_abun_1 = numeric(t_sim),
      total_adult_field_abun_1 = numeric(t_sim),
      total_larvae_feeding_abun_1 = numeric(t_sim),
      total_larvae_pupating_abun_1 = numeric(t_sim),
      average_adult_pat_abun_1 = numeric(t_sim),
      total_adult_pat_abun_2 = numeric(t_sim),
      total_adult_field_abun_2 = numeric(t_sim),
      total_larvae_feeding_abun_2 = numeric(t_sim),
      total_larvae_pupating_abun_2 = numeric(t_sim),
      average_adult_pat_abun_2 = numeric(t_sim),
      total_adult_pat_abun_3 = numeric(t_sim),
      total_adult_field_abun_3 = numeric(t_sim),
      total_larvae_feeding_abun_3 = numeric(t_sim),
      total_larvae_pupating_abun_3 = numeric(t_sim),
      average_adult_pat_abun_3 = numeric(t_sim)
    )
    
    for (i in 1:t_sim) {
      # --- stochastic parameters ------------------------------------------------
      start_biomass <- round(rnorm(1, mean = 1000, sd = 100)) # starting pat biomass
      percent_mortality <- c(1, 0, 0) # 100% adult fly mortality
      # biomass decay rate of pats - full degrad. in six months.
      biomass_decay <- rnorm(10, mean = 5.49, sd = 0.245)
      # spring community structure
      start_abun <- round(rnorm(1, mean = 350, sd = 100) * c(0.5, 0.3, 0.2))
      
      # --- Update pat components and indexes ------------------------------------
      cell_positions <- 1:(size * size) # vector of all cell positions
      # index for cells where there is not a pat
      empty_positions <- which(sapply(pasture, function(cell)
        cell$pat_present == FALSE))
      
      # index for cells where there is a pat
      pat_positions <- which(sapply(pasture, function(cell)
        cell$pat_present == TRUE))
      
      # index for fully degraded pats
      no_biomass <- which(sapply(pasture, function(cell)
        cell$pat_biomass == 0))
      
      # index for IVM treated pats
      treated_positions <-  which(sapply(pasture, function(cell)
        cell$pat_treatment == TRUE))
      
      # index for untreated pats
      untreated_positions <- which(sapply(pasture, function(cell)
        cell$pat_treatment == FALSE))
      
      # pats with no biomass recorded as degraded in pat_present
      for (nob in no_biomass) {
        if (pasture[[nob]]$pat_biomass == 0) {
          pasture[[nob]]$pat_present <- FALSE
        }
        # degraded pats have age reset
        if (pasture[[nob]]$pat_present == FALSE) {
          pasture[[nob]]$pat_age <- NA
        }
      }
      
      # pat age increases by 1 day
      for (age in pat_positions) {
        pasture[[age]]$pat_age <- pasture[[age]]$pat_age + 1
      }
      
      for (pos in empty_positions) {
        # no feeding populations where there isn't a pat
        pasture[[pos]]$larvae_feeding_abun[, ] <- 0
        # no treated status where there isn't a pat
        pasture[[pos]]$pat_treatment <- NA
        
        for (a in 1:guild_no) {
          #' no adult pat populations where there isn't a pat
          pasture[[pos]]$adult_field_abun[a] <- pasture[[pos]]$adult_field_abun[a] + sum(pasture[[pos]]$adult_pat_abun[, a])
          
          pasture[[pos]]$adult_pat_abun[, a] <- 0
        }
      }
      
      # --- update insect populations -------------------------------------------------------------------------
      for (pup in cell_positions) {
        for (h in 1:guild_no) {
          # larvae emerge and join current adults in field
          pasture[[pup]]$adult_field_abun[h] <- pasture[[pup]]$adult_field_abun[h] + pasture[[pup]]$larvae_pupating_abun[pupation_time[h], h]
          
          # emerged larvae removed from pupation df
          pasture[[pup]]$larvae_pupating_abun[pupation_time[h], h] <- 0
        }
        
        # pupating populations age up in all cells
        pasture[[pup]]$larvae_pupating_abun[-1, 1:guild_no] <-
          pasture[[pup]]$larvae_pupating_abun[-nrow(pasture[[pup]]$larvae_pupating_abun), 1:guild_no]
        
        # populations from row 1 removed
        pasture[[pup]]$larvae_pupating_abun[1, 1:guild_no] <- 0
      }
      
      for (fed in pat_positions) {
        # move from feeding to pupating
        for (c in 1:guild_no) {
          pasture[[fed]]$larvae_pupating_abun[1, c] <-
            pasture[[fed]]$larvae_feeding_abun[feeding_time[c], c]
          
          # pupating larvae removed
          pasture[[fed]]$larvae_feeding_abun[feeding_time[c], c] <- 0
          
        }
        
        # feeding populations age up in all cells
        pasture[[fed]]$larvae_feeding_abun[-1, 1:guild_no] <-
          pasture[[fed]]$larvae_feeding_abun[-nrow(pasture[[fed]]$larvae_feeding_abun), 1:guild_no]
        
        #' row 1 populations removed
        pasture[[fed]]$larvae_feeding_abun[1, 1:guild_no] <- 0
      }
      
      for (adu in pat_positions) {
        # adults done in pat move to field
        for (h in 1:guild_no) {
          pasture[[adu]]$adult_field_abun[h] <-
            pasture[[adu]]$adult_field_abun[h] + pasture[[adu]]$adult_pat_abun[(pupation_time[h] + feeding_time[h]), h]
          
          # adults removed from pat if they moved to field
          pasture[[adu]]$adult_pat_abun[(pupation_time[h] + feeding_time[h]), h] <- 0
        }
        # adult pat abun age up in all cells
        pasture[[adu]]$adult_pat_abun[-1, 1:guild_no] <-
          pasture[[adu]]$adult_pat_abun[-nrow(pasture[[adu]]$adult_pat_abun), 1:guild_no]
        
        # row 1 populations removed
        pasture[[adu]]$adult_pat_abun[1, 1:3] <- 0
      }
      
      
      
      # --- Fresh pats generated -------------------------------------------------
      # randomly select fresh pats locations from empty locations
      if (length(empty_positions) == 0) {
        pats_t_I <- empty_positions
      } else {
        pats_t_I <- sample(empty_positions,
                           size = min(pats_per_day, length(empty_positions)),
                           replace = FALSE)
      }
      
      # generate pats and their components in these locations
      for (new in pats_t_I) {
        pasture[[new]]$pat_present <- TRUE
        pasture[[new]]$pat_age <- 0
        pasture[[new]]$pat_biomass <- start_biomass
      }
      # ----------------------------------------------------------------------------
      #' Set starting guild populations
      #' for a controlled amount of days...
      if (i <= start_time) {
        #' flies added to pats generated this time step
        for (new in pats_t_I) {
          pasture[[new]]$adult_pat_abun[1, 1] <- start_abun[1]
        }
        
        #' index for pats up to 5 days old
        beetle_fresh <- which(sapply(pasture, function(cell)
          cell$pat_age <= 5))
        
        #' beetles colonise pats up to 5 days old.
        for (bet in beetle_fresh) {
          pasture[[bet]]$adult_pat_abun[1, 2] <- round(start_abun[2])
          pasture[[bet]]$adult_pat_abun[1, 3] <- round(start_abun[3])
        }
      }
      
      # ----------------------------------------------------------------------------
      #' Determine treatment status of fresh pats
      #' n = number of fresh pats
      #' size = number of trials
      #' prob = chance of being treated (1)
      #' For burn in period (500 days), there is no treatment, then for rest there is.
      if (i < 500) {
        treat_pats <- rep(0, length(pats_t_I))  # All untreated
      } else {
        treat_pats <- rbinom(length(pats_t_I), 1, 1 - tst_level)  # Start binomial draws
      }
      
      #' Assign locations in pats_t_I with corresponding treatment in treat_pats.
      #' Runs the loop for 1 to length(pats_t_I)
      for (tre in seq_along(pats_t_I)) {
        #' gets the position of each pats_t_I
        pos <- pats_t_I[tre]
        #' this returns true for treat_pats = 1, and false for treat_pats = 0.
        pasture[[pos]]$pat_treatment <- treat_pats[tre] == 1
      }
      
      ##############################################################################
      #' Guild adult dispersal
      #' masks for dispersal targets - guild specific
      empty_target <- which(sapply(pasture, function(cell)
        cell$pat_present == FALSE))
      
      fresh_target_beetle <- which(sapply(pasture, function(cell)
        cell$pat_age <= 5))
      fresh_target_fly <- which(sapply(pasture, function(cell)
        cell$pat_age == 0))
      
      old_target_beetle <- which(sapply(pasture, function(cell)
        cell$pat_age > 5))
      old_target_fly <- which(sapply(pasture, function(cell)
        cell$pat_age > 0))
      
      #' guild at a time
      for (d in 1:guild_no) {
        #' empty vector for dispersing populations to each cell
        cumulative_disp <- rep(0, size * size)
        #' index for cells that have dispersing populations
        available_to_disperse <- which(sapply(pasture, function(cell)
          cell$adult_field_abun[d] > 0))
        
        #' for each location...
        for (e in available_to_disperse) {
          #' number of individuals dispersing (for guild d and pat e)
          abundance <- pasture[[e]]$adult_field_abun[d]
          
          #' vector of connectivity values to each destination - check it sums to 1!
          conn_vec <- connectivity_values[e, ]
          conn_vec <- conn_vec / sum(conn_vec) #' should already sum to 1.
          
          #' multinomial distribution to work out population to each destination
          #' n = 1 - each individuals gets 1 destination.
          #' size = abundance - needs to be done for all individuals.
          #' prob = conn_vec - probability that one individuals goes to any cell.
          dispersal_result <- rmultinom(1, size = abundance, prob = conn_vec)
          
          #' accumulate all dispersing individuals for this guild
          cumulative_disp <- cumulative_disp + as.vector(dispersal_result)
          
          #' remove dispersing population from focal cell
          pasture[[e]]$adult_field_abun[d] <- 0
        }
        
        #' dispersal to empty locations for each guild
        for (emp in empty_target) {
          pasture[[emp]]$adult_field_abun[d] <- pasture[[emp]]$adult_field_abun[d] + cumulative_disp[emp]
        }
        
        #' dispersal to fresh locations for guild 1
        if (d == 1) {
          for (fre in fresh_target_fly) {
            pasture[[fre]]$adult_field_abun[d] <- pasture[[fre]]$adult_field_abun[d] + cumulative_disp[fre]
          }
        }
        
        #' dispersal to fresh locations for guilds 2 and 3
        if (d %in% c(2, 3)) {
          for (fre in fresh_target_beetle) {
            pasture[[fre]]$adult_field_abun[d] <- pasture[[fre]]$adult_field_abun[d] + cumulative_disp[fre]
          }
        }
        
        #' dispersal to old locations for guild 1
        if (d == 1) {
          for (old in old_target_fly) {
            pasture[[old]]$adult_field_abun[d] <-  pasture[[old]]$adult_field_abun[d] + cumulative_disp[old]
          }
        }
        
        #' dispersal to old locations for guilds 2 and 3
        if (d %in% c(2, 3)) {
          for (old in old_target_beetle) {
            pasture[[old]]$adult_field_abun[d] <- pasture[[old]]$adult_field_abun[d] + cumulative_disp[old]
          }
        }
      }
      
      #' where dispersal is to a pat, adult population colonises
      #' move fly field populations at pats_t_I into pat
      for (fre in fresh_target_fly) {
        pasture[[fre]]$adult_pat_abun[1, 1] <- pasture[[fre]]$adult_pat_abun[1, 1] +
          pasture[[fre]]$adult_field_abun[1]
        #' remove fly field population
        pasture[[fre]]$adult_field_abun[1] <- 0
      }
      
      #' move beetle field populations at pat_age <5 into pat.
      for (f in 2:3) {
        for (fre in fresh_target_beetle) {
          pasture[[fre]]$adult_pat_abun[1, f] <- pasture[[fre]]$adult_pat_abun[1, f] +
            pasture[[fre]]$adult_field_abun[f]
          #' remove beetle field populations
          pasture[[fre]]$adult_field_abun[f] <- 0
        }
      }
      
      #' index for treated pats
      treated_pats <- which(sapply(pasture, function(cell)
        cell$pat_treatment == TRUE))
      #' treated pats cause mortality for colonising adults
      for (g in 1:guild_no) {
        for (tre in treated_pats) {
          pasture[[tre]]$adult_pat_abun[1, g] <- round(pasture[[tre]]$adult_pat_abun[1, g] * (1 - percent_mortality[g]))
        }
      }
      
      ##############################################################################
      #' Guild oviposition generates feeding population
      #' adult_pat_abun[1,j] in fresh pats goes in , larvae_feeding_abun[1,j] in fresh pats comes out.
      #' Called every time step.
      #' Redo masks just in case
      fresh_target_beetle <- which(sapply(pasture, function(cell)
        cell$pat_age <= 5))
      fresh_target_fly <- which(sapply(pasture, function(cell)
        cell$pat_age == 0))
      
      #' for fresh pats...
      for (fre in fresh_target_fly) {
        #' feeding populations generated based on guild and pat treatment status.
        #' treated pats
        if (pasture[[fre]]$pat_treatment == TRUE) {
          pasture[[fre]]$larvae_feeding_abun[1, 1] <- round(pasture[[fre]]$adult_pat_abun[1, 1] * (fecundity[1, 2] / 2))
        }
        #' untreated pats
        if (pasture[[fre]]$pat_treatment == FALSE) {
          pasture[[fre]]$larvae_feeding_abun[1, 1] <- round(pasture[[fre]]$adult_pat_abun[1, 1] * (fecundity[1, 1] / 2))
        }
      }
      
      #' treated pats
      for (fre in fresh_target_beetle) {
        if (pasture[[fre]]$pat_treatment == TRUE) {
          pasture[[fre]]$larvae_feeding_abun[1, 2] <- round(pasture[[fre]]$adult_pat_abun[1, 2] * (fecundity[2, 2] / 2))
          
          pasture[[fre]]$larvae_feeding_abun[1, 3] <- round(pasture[[fre]]$adult_pat_abun[1, 3] * (fecundity[3, 2] / 2))
        }
        
        #' untreated pats
        if (pasture[[fre]]$pat_treatment == FALSE) {
          pasture[[fre]]$larvae_feeding_abun[1, 2] <- round(pasture[[fre]]$adult_pat_abun[1, 2] * (fecundity[2, 1] / 2))
          
          pasture[[fre]]$larvae_feeding_abun[1, 3] <- round(pasture[[fre]]$adult_pat_abun[1, 3] * (fecundity[3, 1] / 2))
        }
      }
      
      # -----------------------------------------------------------------------------
      #' Dung removal
      for (pos in cell_positions) {
        # determine treatment status for pat
        is_treated <- pasture[[pos]]$pat_treatment &&
          !is.na(pasture[[pos]]$pat_treatment)
        
        # For each guild, get the appropriate dung removal value (treated or untreated)
        # dung_removal is assumed to be a [guild_no x 2] matrix:
        # [,1] = untreated removal per indiv; [,2] = treated removal per indiv
        removal_rates <- if (is_treated)
          dung_removal[, 2]
        else
          dung_removal[, 1]  # length = guild_no
        
        # Get number of feeding individuals per guild
        feeding_inds <- colSums(pasture[[pos]]$larvae_feeding_abun)  # length = guild_no
        
        # Total dung removed per guild (vector), then summed
        total_dung_removed <- sum(removal_rates * feeding_inds)
        
        # Update pat biomass, ensuring it doesn’t go below 0
        pasture[[pos]]$pat_biomass <- max(0, pasture[[pos]]$pat_biomass - total_dung_removed)
        
        # Background decay of dung
        pasture[[pos]]$pat_biomass <- max(0, pasture[[pos]]$pat_biomass - biomass_decay)
      }
      
      ##############################################################################
      #' Update results data frame
      pasture_results[i, ] <- data.frame (
        t_sim = i,
        tst_level = tst_level,
        pats_per_day = pats_per_day,
        total_pat_number = sum(sapply(pasture, function(cell)
          cell$pat_present == TRUE)),
        average_pat_biomass = mean(sapply(pasture, function(cell) {
          if (isTRUE(cell$pat_present)) {
            return(cell$pat_biomass)
          } else {
            return(NA)  # exclude cells where pat_present is not TRUE
          }
        }), na.rm = TRUE),
        total_adult_pat_abun_1 = sum(sapply(pasture, function(cell)
          cell$adult_pat_abun[, 1])),
        total_adult_field_abun_1 = sum(sapply(pasture, function(cell)
          cell$adult_field_abun[1])),
        total_larvae_feeding_abun_1 = sum(sapply(pasture, function(cell)
          cell$larvae_feeding_abun[, 1])),
        total_larvae_pupating_abun_1 = sum(sapply(pasture, function(cell)
          cell$larvae_pupating_abun[, 1])),
        average_adult_pat_abun_1 = sum(sapply(pasture, function(cell)
          sum(cell$adult_pat_abun[, 1]))) / sum(sapply(pasture, function(cell)
            cell$pat_present == TRUE)),
        total_adult_pat_abun_2 = sum(sapply(pasture, function(cell)
          cell$adult_pat_abun[, 2])),
        total_adult_field_abun_2 = sum(sapply(pasture, function(cell)
          cell$adult_field_abun[2])),
        total_larvae_feeding_abun_2 = sum(sapply(pasture, function(cell)
          cell$larvae_feeding_abun[, 2])),
        total_larvae_pupating_abun_2 = sum(sapply(pasture, function(cell)
          cell$larvae_pupating_abun[, 2])),
        average_adult_pat_abun_2 = sum(sapply(pasture, function(cell)
          sum(cell$adult_pat_abun[, 2]))) / sum(sapply(pasture, function(cell)
            cell$pat_present == TRUE)),
        total_adult_pat_abun_3 = sum(sapply(pasture, function(cell)
          cell$adult_pat_abun[, 3])),
        total_adult_field_abun_3 = sum(sapply(pasture, function(cell)
          cell$adult_field_abun[3])),
        total_larvae_feeding_abun_3 = sum(sapply(pasture, function(cell)
          cell$larvae_feeding_abun[, 3])),
        total_larvae_pupating_abun_3 = sum(sapply(pasture, function(cell)
          cell$larvae_pupating_abun[, 3])),
        average_adult_pat_abun_3 = sum(sapply(pasture, function(cell)
          sum(cell$adult_pat_abun[, 3]))) / sum(sapply(pasture, function(cell)
            cell$pat_present == TRUE))
      )
    }
    pasture_results$iteration_id <- tt  # add id column for binding
    pasture_results
  }
stopCluster(cl)
write.csv(treatment_results, "Tables/27-08/test_treatment_results.csv", row.names = FALSE)
