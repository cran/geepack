\name{koch}
\alias{koch}
\non_function{}
\title{ ~~ 1-line description of the data frame ~~ }
\usage{data(koch)}
\description{
The \code{koch} data frame has 288 rows and 4 columns.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{trt}{a numeric vector}
    \item{day}{a numeric vector}
    \item{y}{an ordered factor with levels}
    \item{y}{\code{1} < \code{2} < \code{3}}
    \item{id}{a numeric vector}
  }
}
\details{

}
\source{

}
\references{

}
\examples{
data(koch)
fit <- ordgee(y ~ trt + as.factor(day), id=id, data=koch, corstr="exch")
summary(fit)
}
\keyword{datasets}
