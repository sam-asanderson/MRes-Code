#Fly LV model - single breeding/dispersal event
# Parameters
N0_f <- 0          # Initial population size
r_f <- 0.099          # Growth rate - Vale and Grant (2002) fast breeding species
K_f <- 559         # Carrying capacity (all fly species - total flies collected/total number of pats)
T_f <- 60           # Number of time steps
I <- 100          #Immigration 
# Initialize population array
N_f <- numeric(T_f + 1)
N_f[1] <- N0_f

# Run the model
for (t in 1:1) {
  N_f[t + 1] <- N_f[t] + I
}
for (t in 2:60) {
  N_f[t] + 1] <- N_f[t] + r_f * N_f[t] * (1 - N_f[t] / K_f)
}

#plot results
plot(0:T_f, N_f, type = "o", pch = 19, ylim = c(0, 600), col = "blue", xlab = "Time step", ylab = "Population size",
     main = "Single-Species Discrete Lotka-Volterra Model for Flies")
