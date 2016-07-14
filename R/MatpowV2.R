Matpow <- function(M, numer, denom = 1, digits = 2){
  # Raises any diagonalizable Matrix to any power (even complex).
  # Eigenvectors can be complex
  # Capabilities for non-diagonalizable matrices:
  #   Exponentiation (integer powers only)
  #   Square root
  #
  # Args:
  #   M: a square Matrix
  #   numer: numerator of (rational) exponent. Can be a decimal.
  #   denom: denominator of rational exponent (1 by default)
  #
  # Returns:
  #   The solution to the exponentiation operation supplied
  #   Allows to take powers and roots
  #   Allows to compute the inverse matrix (if invertible)
  #   Method based on spectral decomposition
  #   Returns a real-valued root whenever possible
  #   Returns the (principal) complex root if that is the only root
  frac <- (numer / denom)
  e.val <- eigen(M)$values
  n <- length(e.val)
  V <- eigen(M)$vectors
  D <- diag(eigen(M)$values)
  dexp <- numeric(n)
  if (rankMatrix(V) != n){
    if(frac == 0.5){
      sqrtm(M)
    }else{
      if(numer %% 1 == 0){
        M %^% numer
      }else{
        print("Sorry, requested operation requires a diagonalizable matrix.")
      }
    }
  }else{
    Im.zap <- function(x) {
      if (all(Im(z <- zapsmall(x)) == 0)) as.numeric(z) else x
    }
    if ((rankMatrix(M) != n) && frac == -1){
      stop("Sorry, requested operation requires a non-singular matrix.")
    }else{
    for (i in 1:n){
      dexp[i] <- explus(e.val[i], numer, denom)
    }
    Dexp <- diag(dexp)
    Mexp <- V %*% Dexp %*% solve(V)
    solution <- round(matrix(Im.zap(Mexp), ncol = n), digits)
    return(solution)
    }
  }
}
