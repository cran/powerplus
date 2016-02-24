\name{Matpow}
\alias{Matpow}
\title{Matrix Power}
\description{
Raises any diagonalizable Matrix to any real power.
Capabilities for non-diagonalizable matrices:
Exponentiation (integer powers only), Square root
}
\usage{
Matpow(M, numer, denom = 1)
}
\arguments{
  \item{M}{a square Matrix}
  \item{numer}{numerator of (rational) exponent. Can be a decimal number.}
  \item{denom}{denominator of rational exponent (1 by default)}
}
\details{
Method based on spectral decomposition.
Allows to take powers and roots.
Allows to compute the inverse matrix (if invertible).
}
\value{
The solution to the exponentiation operation supplied.
Computes a real-valued root whenever possible. Otherwise, the principal complex root.
}
\references{
For more on spectral decomposition (also known as eigendecomposition), visit \url{http://mathworld.wolfram.com/EigenDecomposition.html}
}
\author{Albert Dorador}
\examples{
A <- matrix(1:4, ncol = 2)
Matpow(A, 3)
Matpow(A, 0.5)
Matpow(A, 0.2)
Matpow(A, 1, 5)
Matpow(A, 2, 4)
Matpow(A, -1)

C <- matrix(sample(1:1000, 81), ncol = 9)
Matpow(C, 2)
Matpow(C, 0.5)
Matpow(C, 4, 5)
Matpow(C, pi)
Matpow(C, 0.73)
Matpow(C, -1)

M <- matrix(c(1, 0, 1, 1), 2, 2)  # A non-diagonalizable matrix
Matpow(M, 3)
Matpow(M, 0.5)
Matpow(M, 4, 8)
}
