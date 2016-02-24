\name{explus}
\alias{explus}
\title{Scalar Exponentiation}
\description{Raises any base (real or complex) to any real power.}
\usage{
explus(a, numer, denom = 1)
}
\arguments{
  \item{a}{any base (real or complex)}
  \item{numer}{numerator of (rational) exponent. Can be a decimal number.}
  \item{denom}{denominator of rational exponent (1 by default)}
}
\value{
The solution to the exponentiation operation supplied.
Returns a real-valued root whenever possible. Otherwise, the principal complex root.
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
}
