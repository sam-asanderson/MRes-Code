rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Initialize the simulation parameters
n_rows <- 10  # number of rows in the grid
n_cols <- 10  # number of columns in the grid
day_counter <- 0
num_pats_to_add <- 10

# Create an empty grid with patch properties
grid <- expand.grid(row = 1:n_rows, col = 1:n_cols)
grid$resources <- 0
grid$pat_status <- "none"
grid$pat_lifetime <- 0
grid$pcolor <- "green"

# Function to increment the day counter
increment_day <- function() {
  day_counter <<- day_counter + 1
}

# Function to add new pats to the grid
add_new_pats <- function() {
  available_patches <- grid %>% filter(pat_status == "none")
  
  if (nrow(available_patches) > 0) {
    for (i in 1:num_pats_to_add) {
      if (nrow(available_patches) == 0) break  # Stop if no available patches
      
      target_patch_index <- sample(1:nrow(available_patches), 1)
      target_patch <- available_patches[target_patch_index, ]
      
      if (runif(1) < 0.5) {
        grid[grid$row == target_patch$row & grid$col == target_patch$col, ] <- target_patch %>%
          mutate(pat_status = "treated", pat_lifetime = 100, resources = 500, pcolor = "red")
      } else {
        grid[grid$row == target_patch$row & grid$col == target_patch$col, ] <- target_patch %>%
          mutate(pat_status = "untreated", pat_lifetime = 50, resources = 500, pcolor = "brown")
      }
    }
  }
}

# Function to update patches
update_pats <- function() {
  for (i in 1:nrow(grid)) {
    if (grid$pat_status[i] != "none") {
      grid$pat_lifetime[i] <- grid$pat_lifetime[i] - 1
      if (grid$pat_lifetime[i] <= 0) {
        grid$pat_status[i] <- "none"
        grid$resources[i] <- 0
        grid$pcolor[i] <- "green"
      }
    }
  }
}

# Set up the initial state
setup <- function() {
  grid$pcolor <- "green"
  grid$pat_status <- "none"
  grid$pat_lifetime <- 0
  grid$resources <- 0
}

# Run one step of the simulation
go <- function() {
  increment_day()
  add_new_pats()
  update_pats()
}

# Function to plot the grid
plot_grid <- function() {
  ggplot(grid, aes(x = col, y = row, fill = pcolor)) +
    geom_tile() +
    scale_fill_manual(values = c("green" = "green", "red" = "red", "brown" = "brown")) +
    theme_minimal() +
    coord_fixed() +
    labs(title = paste("Day:", day_counter))
}

# Set up initial grid state
setup()

# Simulate for 30 days (for example)
for (day in 1:30) {
  go()  # Run the simulation step
  plot_grid()  # Plot the grid at each step
}
