\name{respdis}
\alias{respdis}
\non_function{}
\title{Clustered Ordinal Respiratory Disorder}
\usage{data(respdis)}
\description{
The \code{respdis} data frame has 444 rows and 3 columns. The study
described in Miller et. al. (1993) is a randomized clinical trial of a
new treatment of respiratory disorder.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{a factor  for subject id with levels 0 - 110}
    \item{resp}{an ordered factor for the response with levels}
    \item{resp}{\code{1} < \code{2} < \code{3}, 1 = poor, 2 = good, and
      3 = excellent}
    \item{trt}{a factor for treatment with levels}
    \item{trt}{\code{-1} }
    \item{trt}{\code{1} }
  }
}
\references{
Miller, M.E., David, C.S., and Landis, R.J. (1993) The analysis of
longitudinal polytomous data: Generalized estimating equation and
connections with weighted least squares, \emph{Biometrics} \bold{49}: 1033-1048.
}
\examples{
data(respdis)
fit <- ordgee(resp ~ trt, id=id, data=respdis, int.const=FALSE)
summary(fit)
}
\keyword{datasets}
