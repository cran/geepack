\name{seizure}
\alias{seizure}
\non_function{}
\title{Epiliptic Seizures}
\usage{data(seizure)}
\description{
The \code{seizure} data frame has 295 rows and 5 columns. The dataset
has the number of epiliptic seizures in each of four two-week intervals,
for treatment and control groups with a total of 59 individuals. Visit 0
is the baseline measured in an interval of 8 weeks. 
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{y}{a numeric vector of the number of epiliptic seizures}
    \item{trt}{an indicator of treatment}
    \item{id}{a numeric vector of subject id}
    \item{visit}{a numeric vector of visit number}
    \item{age}{a numeric vector of subject age}
  }
}
\references{
  Diggle, P.J., Liang, K.Y., and Zeger, S.L. (1994) Analysis of Longitudinal
  Data, Clarendon Press.
}
\examples{
data(seizure)
seizure$t <- ifelse(seizure$visit == 0, 8, 2)
seizure$x <- ifelse(seizure$visit == 0, 0, 1)
mm <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data=seizure, corstr="exch", family=poisson)
summary(mm)
}
\keyword{datasets}
