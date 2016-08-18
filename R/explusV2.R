#' @title Scalar Exponentiation
#' @description Raises any base (real or complex) to any power (even complex).
#' @details Method based on numerical treatment of complex exponents using Euler Formula and angles
#' measured in radians, which is the default method built-in in R.
#' @references
#' For more on Euler Formula, visit
#' \url{http://mathworld.wolfram.com/EulerFormula.html}
#'
#' For more on complex exponents, visit
#' \url{http://mathworld.wolfram.com/ComplexExponentiation.html}
#' @param a any base (real or complex).
#' @param numer numerator of exponent. Can be a decimal or complex number.
#' @param denom denominator of exponent (1 by default). Can be a decimal or complex number.
#' @return The solution to the exponentiation operation supplied. Returns a real-valued root whenever possible.
#' Otherwise, the principal complex root.
#' @author Albert Dorador
#' @export
#' @import phonTools
#' @importFrom Matrix rankMatrix
#' @import expm
#' @import MASS
#' @examples
#' explus(-3, 4, 2)
#' explus(-3, 2, 4)
#' explus(-3, 2, 3)
#' explus(-3, 5, 3)
#' explus(-3, 5, 2)
#' explus(-3, -2, 4)
#' explus(0-0.5773503i, 2)
#' explus(-0.4, pi)
#' explus(-0.37, 0.2)
#' explus(-0.37, 1, 5)
#' explus(5, 7i)
#' explus(2+3i, 1+2i)
#' explus(2+3i, 1+2i, -4+1i)
#' explus(2+3i, 1+2i, 8)
#'

explus <- function(a, numer, denom = 1){
  if (denom == 0)
    stop("Dividing by 0 breaks mathematics")
  if (a == 0){
    if (numer == 0){
      return(1) #by R convention
    }else{
      return(0)
    }
  }
  if (is.complex(numer / denom)){
    a ** (numer / denom)
  }else{
    if ((numer %% 1 == 0) && (denom %% 1 == 0)){ # Case 0: numerator is integer, denominator is integer
      vec <- reduce.fraction (c(numer, denom))
      numer <- vec[1]
      denom <- vec[2]
      if (is.complex(a)){
        a ** (numer / denom)
      }else{
        if (a < 0){
          if (numer %% 2 == 0){
            (a ** numer) ** (1 / denom)
          }else{
            if (denom %% 2 == 0){
              (as.complex(a)) ** (numer / denom)
            }else{
              Nthroot <- function(x, n){ #only for n odd
                sign(x) * abs(x) ** (1 / n)
              }
              (Nthroot(a, denom)) ** numer
            }
          }
        }else{
          a ** (numer / denom)
        }
      }
    }else{ # There exist 3 possible mutually exclusive cases
      if ((numer %% 1 != 0) && (denom %% 1 == 0)){ # Case 1: numerator is decimal, denominator is integer
        denom_of_numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[2]
        numer_of_numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[1]
        numer <- numer_of_numer
        denom <- denom * denom_of_numer
      }
      if ((numer %% 1 == 0) && (denom %% 1 != 0)){ # Case 2: numerator is integer, denominator is decimal
        denom_of_denom <- as.numeric(unlist(strsplit(attributes(fractions(denom))$fracs,split="/")))[2]
        numer_of_denom <- as.numeric(unlist(strsplit(attributes(fractions(denom))$fracs,split="/")))[1]
        numer <- numer * denom_of_denom
        denom <- numer_of_denom
      }
      if ((numer %% 1 != 0) && (denom %% 1 != 0)){ # Case 3: numerator is decimal, denominator is decimal
        denom_of_numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[2]
        numer_of_numer <- as.numeric(unlist(strsplit(attributes(fractions(numer))$fracs,split="/")))[1]
        denom_of_denom <- as.numeric(unlist(strsplit(attributes(fractions(denom))$fracs,split="/")))[2]
        numer_of_denom <- as.numeric(unlist(strsplit(attributes(fractions(denom))$fracs,split="/")))[1]
        numer <- numer_of_numer * denom_of_denom
        denom <- denom_of_numer * numer_of_denom
      }
      if (is.complex(a)){
        a ** (numer / denom)
      }else{
        if (a < 0){
          if (numer %% 2 == 0){
            (a ** numer) ** (1 / denom)
          }else{
            if (denom %% 2 == 0){
              (as.complex(a)) ** (numer / denom)
            }else{
              Nthroot <- function(x, n){ #only for n odd
                sign(x) * abs(x) ** (1 / n)
              }
              (Nthroot(a, denom)) ** numer
            }
          }
        }else{
          a ** (numer / denom)
        }
      }
    }
  }
}
