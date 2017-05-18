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
#' @param a any base (real or complex). It can be a vector.
#' @param numer numerator of exponent. Can be a decimal or complex number.
#' @param denom denominator of exponent (1 by default). Can be a decimal or complex number.
#' @param n.cycles The maximum number of steps to be used in the continued fraction approximation process done internally.
#' Default is 10. Increasing this default may increase precision, but it may also produce unexpected results (even a runtime error).
#' @param tol a tolerance, \eqn{10^{-12}}{10^-12} by default. Prevents possible numerical problems.
#' Can be set to 0 if desired.
#' @return The solution to the exponentiation operation supplied. Returns a real-valued root whenever possible.
#' Otherwise, the principal complex root.
#' @author Albert Dorador
#' @export
#' @import complexplus
#' @import phonTools
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

explus <- function(a, numer, denom = 1, n.cycles = 10, tol = 1e-12){
  numer <- Imzap(numer, tol) #just in case the user tries e.g. 2+0i
  denom <- Imzap(denom, tol) #same idea
  if (denom == 0)
    stop("Dividing by 0 breaks mathematics")
  lvec = length(a)
  solution = numeric(lvec)
  for (i in 1:lvec){
    if (a[i] == 0){
      if (numer == 0){
        solution[i] = 1 #by R convention
      }else{
        solution[i] = 0
      }
    }
    if (is.complex(numer / denom)){
      solution[i] = a[i] ** (numer / denom)
    }else{
      if ((numer %% 1 == 0) && (denom %% 1 == 0)){ # Case 0: numerator is integer, denominator is integer
        vec <- reduce.fraction (c(numer, denom))
        numer <- vec[1]
        denom <- vec[2]
        if (is.complex(a[i])){
          solution[i] = a[i] ** (numer / denom)
        }else{
          if (a[i] < 0){
            if (numer %% 2 == 0){
              solution[i] =  (a[i] ** numer) ** (1 / denom)
            }else{
              if (denom %% 2 == 0){
                solution[i] = (as.complex(a[i])) ** (numer / denom)
              }else{
                Nthroot <- function(x, n){ #only for n odd
                  sign(x) * abs(x) ** (1 / n)
                }
                solution[i] = (Nthroot(a[i], denom)) ** numer
              }
            }
          }else{
            solution[i] = a[i] ** (numer / denom)
          }
        }
      }else{ # There exist 3 possible cases; making sure numer & denom of one case do not pass to the next!
        if ((numer %% 1 != 0) && (denom %% 1 == 0)){ # Case 1: numerator is decimal, denominator is integer
          denom_of_numer <- as.double(unlist(strsplit(attributes(fractions(numer,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[2]
          numer_of_numer <- as.double(unlist(strsplit(attributes(fractions(numer,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[1]
          numer <- numer_of_numer
          denom <- denom * denom_of_numer
        }
        else if ((numer %% 1 == 0) && (denom %% 1 != 0)){ # Case 2: numerator is integer, denominator is decimal
          denom_of_denom <- as.double(unlist(strsplit(attributes(fractions(denom,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[2]
          numer_of_denom <- as.double(unlist(strsplit(attributes(fractions(denom,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[1]
          numer <- numer * denom_of_denom
          denom <- numer_of_denom
        }
        else { # Case 3: numerator is decimal, denominator is decimal
          denom_of_numer <- as.double(unlist(strsplit(attributes(fractions(numer,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[2]
          numer_of_numer <- as.double(unlist(strsplit(attributes(fractions(numer,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[1]
          denom_of_denom <- as.double(unlist(strsplit(attributes(fractions(denom,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[2]
          numer_of_denom <- as.double(unlist(strsplit(attributes(fractions(denom,cycles=n.cycles,max.denominator=Inf))$fracs,split="/")))[1]
          numer <- numer_of_numer * denom_of_denom
          denom <- denom_of_numer * numer_of_denom
        }
        if (is.complex(a[i])){
          solution[i] = a[i] ** (numer / denom)
        }else{
          if (a[i] < 0){
            if (numer %% 2 == 0){
              solution[i] = (a[i] ** numer) ** (1 / denom)
            }else{
              if (denom %% 2 == 0){
                solution[i] = (as.complex(a[i])) ** (numer / denom)
              }else{
                Nthroot <- function(x, n){ #only for n odd
                  sign(x) * abs(x) ** (1 / n)
                }
                solution[i] = (Nthroot(a[i], denom)) ** numer
              }
            }
          }else{
            solution[i] = a[i] ** (numer / denom)
          }
        }
      }
    }
  }
  return(Imzap(solution, tol))
}
