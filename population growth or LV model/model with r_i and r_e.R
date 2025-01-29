#Flies
N=rep(0,60); N[1]=10;
r_b <- 0.219
r_d <- 0.12
r_i <- (0.4 * exp(-0.4 * N[i])) #taken from Hanski, 1980.#
r_e <- (0.15 * (1 - exp(-0.04 * N[i]))) #"
K_f <- 559 # adding carrying capacity doesn't matter here - (1 - N[i] / K_f).

# Iterate to compute N
for (i in 1:59) {
  N[i+1]=N[i] + (r_b * N[i] + (0.4 * exp(-0.4 * N[i])) * N[i] - r_d * N[i] - (0.15 * (1 - exp(-0.04 * N[i]))) * N[i])
}

plot(
  0:(length(N)-1), N, 
  type = "o",  # Line with points
  col = "blue",  # Line color
  pch = 16,  # Point style
  xlab = "Time Step", 
  ylab = "Population Size",
  main = "Fly population Dynamics Over Time", 
  las = 1,  # Rotate y-axis labels for readability
  bty = "l"  # Box around plot
)
grid()  # Add a grid
 