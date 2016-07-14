\name{explus}
\alias{explus}
\title{Scalar Exponentiation}
\description{Raises any base (real or complex) to any power (even complex).}
\usage{
explus(a, numer, denom = 1)
}
\arguments{
  \item{a}{any base (real or complex)}
  \item{numer}{numerator of exponent. Can be a decimal or complex number.}
  \item{denom}{denominator of exponent (1 by default). Can be a decimal or complex number.}
}
\details{
Method based on numerical treatment of complex exponents using Euler's Formula and angles measured in radians, which is the default method built-in in R.
}
\value{
The solution to the exponentiation operation supplied.
Returns a real-valued root whenever possible. Otherwise, the principal complex root.
}
\references{
For more on Euler's Formula, visit \url{http://mathworld.wolfram.com/EulerFormula.html}

For more on complex exponents, visit \url{http://mathworld.wolfram.com/ComplexExponentiation.html}
}
\author{Albert Dorador}
\examples{
explus(-3, 4, 2)
explus(-3, 2, 4)
explus(-3, 2, 3)
explus(-3, 5, 3)
explus(-3, 5, 2)
explus(-3, -2, 4)
explus(0-0.5773503i, 2)
explus(-0.4, pi)
explus(-0.37, 0.2)
explus(-0.37, 1, 5)
explus(5, 7i)
explus(2+3i, 1+2i)
explus(2+3i, 1+2i, -4+1i)
explus(2+3i, 1+2i, 8)
}
