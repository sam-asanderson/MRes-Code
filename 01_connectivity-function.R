#' Generate Patch Connectivities
#'
#' Construct a matrix of connectivities for a given landscape motif and number of patches
#'
#' @param size Numeric. Number of patches.
#' @param rate Numeric. If `fall.off = "exponential"`, the rate parameter of the exponential function.
#' @returns A symmetrical square matrix of pairwise patch connectivities. Column is focal cell number, row is target cell number. 
#'
#' @examples
#' P <- generateP(size = 25, norm = TRUE, rate = 1)
#'  

GenerateConn <- function(size, norm = TRUE, rate) {
  # creates numbered size * size matrix.]
  mat <- matrix(1:(size * size),
                ncol = size,
                nrow = size,
                byrow = FALSE)
  
  # Gives all cells an ID 1:size.
  pos <- arrayInd(1:(size * size), .dim = c(size, size))
  
  # Creates expanded matrix - all cells * all cells.
  dist_mat <- matrix(0, nrow = (size * size), ncol = (size * size))
  
  # Compute Chebyshev (8 directions) distances between all pairs - abs is absolute value.
  for (x in 1:(size * size)) {
    for (y in x:(size * size)) {
      d <- max(abs(pos[x, 1] - pos[y, 1]), abs(pos[x, 2] - pos[y, 2]))
      
      # Fill in dist_mat with distance values - symmetrical connectivity.
      dist_mat[x, y] <- d
      dist_mat[y, x] <- d
    }
  }
  
  # Generate exponential fall-off connectivity values - is normalization here needed?
  d_vals <- 1:(size - 1)
  probs <- (1 - pexp(d_vals, rate))
  probs <- probs / sum(probs)
  
  # Apply connectivity values to all cells.
  for (r in 1:(size - 1)) {
    dist_mat[dist_mat == r] <- probs[r]
  }
  
  # Make the matrix a row-stochastic matrix - now shows dispersal probability.
  if (norm) {
    dist_mat <- sweep(dist_mat, 1, rowSums(dist_mat), FUN = "/")
  }
  
  return(dist_mat)
}


