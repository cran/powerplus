#' @title Matrix Power
#' @description Raises a valid Matrix to any power (even complex).
#' Valid matrices are square matrices that are diagonalizable or whose real eigenvalues are positive.
#' @details If the matrix is diagonalizable, the method used is based on spectral decomposition;
#' if the matrix is not diagonalizable, the method used is based on matrix exponentials and logarithms,
#' calling functions \code{matexp} and \code{matlog}, both from package \pkg{complexplus}.
#' The particular method used to compute the matrix exponential and logarithm may be chosen from the options available
#' in functions \code{expm} and \code{logm} respectively, both from package \pkg{expm}.
#' Note that \code{Matpow}, by extension, allows one to compute roots and the matrix inverse (if invertible).
#' @references
#' For more on spectral decomposition (also known as eigendecomposition), visit
#' \url{http://mathworld.wolfram.com/EigenDecomposition.html}
#' @param M a square matrix
#' @param numer numerator of exponent. Can be a decimal or complex number.
#' @param denom denominator of exponent (1 by default). Can be a decimal or complex number.
#' @param expmethod method chosen to compute the matrix exponential if matrix is known to be
#' non-diagonalizable. The default method is the same as in function \code{expm} from package \pkg{expm}.
#' @param logmethod method chosen to compute the matrix logarithm if matrix is known to be
#' non-diagonalizable. The default method is the same as in function \code{logm} from package \pkg{expm}.
#' @param tol a tolerance, \eqn{10^{-12}}{10^-12} by default. Prevents possible numerical problems.
#' Can be set to 0 if desired.
#' @return The solution to the exponentiation operation supplied.
#' For diagonalizable matrices, \code{Matpow} returns a real-valued root whenever possible
#' (otherwise, the principal complex root).
#' @author Albert Dorador
#' @export
#' @import complexplus
#' @importFrom Matrix rankMatrix
#' @import expm
#' @seealso
#' \code{\link[complexplus]{matexp}}
#' \code{\link[complexplus]{matlog}}
#' \code{\link[expm]{expm}}
#' \code{\link[expm]{logm}}
#' @examples
#' A <- matrix(1:4, ncol = 2)
#'
#' Matpow(A, 3)
#' Matpow(A, 0.5)
#' Matpow(A, 0.2)
#' Matpow(A, 1, 5)
#' Matpow(A, 2, 4, expmethod = "Pade", logmethod = "Eigen") #inocuous, as A is diagonalizable
#' Matpow(A, -1)
#' Matpow(A, 2+5i)
#' Matpow(A, 3i)
#' Matpow(A, 1+2i)
#' Matpow(A, 3i, 2+7i)
#'
#' B <- matrix(sample(1:100, 81), ncol = 9)
#' Matpow(B, 2)
#' Matpow(B, 0.5)
#' Matpow(B, 7+2i)
#' Matpow(B, 4i, 1+3i)
#'
#' C <- matrix(c(1, 0, 1, 1), ncol = 2) # A non-diagonalizable matrix
#' Matpow(C, 3)
#' Matpow(C, 0.5)
#' Matpow(C, 4, 8, expmethod = "Taylor", logmethod = "Eigen")
#' Matpow(C, 0.5*pi)
#' Matpow(C, 0.24)
#' Matpow(C, -2)
#' Matpow(C, 3+5i)
#' Matpow(C, 2i, 1+9i)
#'

Matpow <- function(M, numer, denom = 1, expmethod = "Higham08.b", logmethod = "Higham08", tol = 1e-12){
  if (denom == 0)
    stop("Dividing by 0 breaks mathematics")
  d <- dim(M)
  n <- d[1]
  p <- d[2]
  if (n != p)
    stop("Supplied matrix is not square")
  numer <- Imzap(numer, tol) #just in case the user tries e.g. 2+0i
  denom <- Imzap(denom, tol) #same idea
  power <- (numer / denom)
  eigen <- eigen(M)
  e.val <- eigen$values
  V <- eigen$vectors
  D <- diag(e.val)
  dexp <- numeric(n)
  
  ## Always attempt sqrtm; conjecture: sqrtm is fully optimized ##
  if(power == 0.5){
    solution <- matrix(Imzap(sqrtm(M), tol), ncol = n)
    return(solution)
  }else{
    
    ## ONLY for NON-Complex powers ##
    if (!is.complex(power)){  #if power is not complex (%% and < are undefined for complex)
      if ((rankMatrix(M) != n) && (power < 0)){ #check if she wants to invert a non-invertible matrix
        stop("Requested operation requires a non-singular matrix.")
      }else{
        if((power %% 1 == 0) && (power >= 0)){ #because power is NOT complex, try to use %^%
          solution <- matrix(Imzap(M %^% power, tol), ncol = n)
          return(solution)
        }else{ #if power is not complex but is NOT a non-neg integer, e.g. power is -8 or 9.3
          if (rankMatrix(V) != n){  # ONLY if Matrix is NON-Diagonalizable #
            safe <- TRUE
            i <- 1
            while (safe && i <= n){
              if(class(Imzap(e.val[i], tol)) != "complex"){ #i.e. if ith eigenvalue is real
                if(Imzap(e.val[i], tol) <= 0){ #0 would make the matrix singular, and negative is likely to cause trouble
                  safe <- FALSE
                }
              }
              i <- i + 1
            }
            if (safe){ #Go ahead: taking the log is safe
              solution <- matrix(Imzap(matexp(power * matlog(M, logmethod), expmethod), tol), ncol = n)
              return(solution) #I'm aware of slight inefficiency in repeating some code in matlog
            }else{
              print("Sorry, requested operation requires a diagonalizable matrix.")
            }
          }else{  # if matrix IS Diagonalizable #
            dexp <- vapply(e.val, function(a, b, c = 1) explus(a, numer, denom, tol=tol), complex(1))
            Dexp <- diag(dexp)
            Mexp <- V %*% Dexp %*% solve(V)
            solution <- matrix(Imzap(Mexp, tol), ncol = n)
            return(solution)
          }
        }
      }
      
      ## ONLY for Complex powers ##
    }else{
      if (rankMatrix(V) != n){  # ONLY if Matrix is NON-Diagonalizable #
        safe <- TRUE
        i <- 1
        while (safe && i <= n){
          if(class(Imzap(e.val[i], tol)) != "complex"){ #i.e. if ith eigenvalue is real
            if(Imzap(e.val[i], tol) <= 0){ #0 would make the matrix singular, and negative is likely to cause trouble
              safe <- FALSE
            }
          }
          i <- i + 1
        }
        if (safe){ #Go ahead: taking the log is safe
          solution <- matrix(Imzap(matexp(power * matlog(M, logmethod), expmethod), tol), ncol = n)
          return(solution) #I'm aware of slight inefficiency in repeating some code in matlog
        }else{
          print("Sorry, requested operation requires a diagonalizable matrix.")
        }
      }else{  # if matrix IS Diagonalizable #
        dexp <- vapply(e.val, function(a, b, c = 1) explus(a, numer, denom, tol=tol), complex(1))
        Dexp <- diag(dexp)
        Mexp <- V %*% Dexp %*% solve(V)
        solution <- matrix(Imzap(Mexp, tol), ncol = n)
        return(solution)
      }
    }
  }
}
