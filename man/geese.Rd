\name{geese}
\alias{geese}
\alias{geese.fit}
\alias{print.geese}
\alias{summary.geese}
\alias{print.summary.geese}
\title{Function to solve a Generalize Estimating Equation Model}
\description{
  Produces an object of class `geese' which is a Generalized Estimating
  Equation fit of the data.
}
\usage{
geese(formula = formula(data), sformula = ~1, id, waves = NULL,
      data = parent.frame(), subset = NULL, na.action = na.omit,
      contrasts = NULL, weights = NULL, zcor = NULL, corp = NULL,
      control = geese.control(...), b = NA, alpha = NA, gm = NA,
      family = gaussian(), mean.link = NULL, variance = NULL,
      cor.link = "identity", sca.link = "identity", link.same = TRUE,
      scale.fix = F, scale.value = 1, corstr = "independence", ...)

geese.fit(x, y, id, offset = rep(0, N), soffset = rep(0, N),
          weights = rep(1,N), waves = NULL, zsca = matrix(1, N, 1),
          zcor = NULL, corp = NULL, control = geese.control(...),
          b = NA, alpha = NA, gm = NA, family = gaussian(),
          mean.link = NULL, variance = NULL, cor.link = "identity",
          sca.link = "identity", link.same = TRUE, scale.fix = F,
          scale.value = 1, corstr = "independence", ...)
}
\arguments{
  \item{formula}{a formula expression as for \code{glm}, of the form
    \code{response ~ predictors}. See the documentation of lm and
    formula for details. As for glm, this specifies the linear predictor
    for modeling the mean. A term of the form \code{offset(expression)}
    is allowed.
  }
  \item{sformula}{a formula expression of the form \code{ ~ predictor},
    the response being ignored. This specifies the linear predictor for
    modeling the dispersion. A term of the form
    \code{offset(expression)} is allowed.
  }
  \item{id}{a vector which identifies the clusters.  The length of `id'
    should be the same as the number of observations.  Data are assumed
    to be sorted so that observations on a cluster are contiguous rows
    for all entities in the formula.
  }
  \item{waves}{an integer vector which identifies components in
    clusters. The length of \code{waves} should be the same as the
    number of observation. components with the same \code{waves} value
    will have the same link functions. 
  }
  \item{data}{
    an optional data frame in which to interpret the variables occurring
    in the \code{formula}, along with the \code{id} and \code{n} variables.
  }
  \item{subset}{
    expression saying which subset of the rows of the data should be used
    in the fit.  This can be a logical vector (which is replicated to have
    length equal to the number of observations), or a numeric vector
    indicating which observation numbers are to be included, or a
    character vector of the row names to be included.
    All observations are included by default.
  }
  \item{na.action}{
    a function to filter missing data.  For \code{gee} only \code{na.omit}
    should be used here.
  }
  \item{contrasts}{
    a list giving contrasts for some or all of the factors appearing
    in the model formula.  The elements of the list should have the
    same name as the variable and should be either a contrast matrix
    (specifically, any full-rank matrix with as many rows as there are
    levels in the factor), or else a function to compute such a matrix
    given the number of levels.
  }
  \item{weights}{an optional vector of weights to be used
    in the fitting process. The length of \code{weights} should be the
    same as the number of observations.
  }
  \item{zcor}{a design matrix for correlation parameters.}
  \item{corp}{known parameters such as coordinates used for correlation
    coefficients.
  }
  \item{control}{a list of iteration and algorithmic constants. See
    \code{\link{geese.control}} for their names and default
    values. These can also be set as arguments to \code{geese} itself.
  }
  \item{b}{an initial estimate for the mean parameters.}
  \item{alpha}{an initial estimate for the correlation parameters.}
  \item{gm}{an initial estimate for the scale parameters.}
  \item{family}{a description of the error distribution and link
    function to be used in the model, as for \code{\link{glm}}.
  }
  \item{mean.link}{a character string specifying the link function for
    the means. The following are allowed:
    \code{"identity"}, \code{"logit"}, \code{"probit"},
    \code{"cloglog"}, \code{"log"}, and \code{"inverse"}.
    The default value is determined from family.
  }
  \item{variance}{a character string specifying the variance function
    in terms of the mean. The following are allowed:
    \code{"gaussian"}, \code{"binomial"}, \code{"poisson"}, and
    \code{"gamma"}. The default value is determined from family.
  }
  \item{cor.link}{a character string specifying the link function for
    the correlation coefficients. The following are allowed:
    \code{"identity"}, and \code{"fisherz"}.
  }
  \item{sca.link}{a character string specifying the link function for
    the scales. The following are allowed:
    \code{"identity"}, and \code{"log"}.
  }
  \item{link.same}{a logical indicating if all the components in a
    cluster should use the same link.
  }
  \item{scale.fix}{
    a logical variable; if true, the scale parameter is fixed at
    the value of \code{scale.value}.
  }
  \item{scale.value}{
    numeric variable giving the value to which the scale parameter
    should be fixed; used only if \code{scale.fix == TRUE}.
  }
  \item{corstr}{
    a character string specifying the correlation structure.
    The following are permitted:
    \code{"independence"},
    \code{"exchangeable"},
    \code{"ar1"},
    \code{"unstructured"} and
    \code{"userdefined"}
  }
  \item{x, y}{
    \code{x} is a design matrix of dimension \code{n * p}, and \code{y}
    is a vector of observations of length \code{n}.
  }
  \item{offset, soffset}{
    vector of offset for the mean and for the scale, respectively.
  }
  \item{zsca}{
    a design matrix of dimension \code{n * r} for the scales.
  }
  \item{\dots}{further arguments passed to or from other methods.}
}
\value{
  An object of class \code{"geese"} representing the fit.
}
\references{
  Liang, K.Y. and Zeger, S.L. (1986)
  Longitudinal data analysis using generalized linear models.
  \emph{Biometrika}, \bold{73} 13--22.

  Prentice, R.L. and Zhao, L.P. (1991).
  Estimating equations for parameters in means and covariances of
  multivariate discrete and continuous responses.
  \emph{Biometrics}, \bold{47} 825--839.
}
\author{Jun Yan \email{jyan@stat.wisc.edu}}
\seealso{
  \code{\link{glm}}, \code{\link{lm}}, \code{\link{ordgee}}.
}
\examples{
data(seizure)
seizure$t <- ifelse(seizure$visit == 0, 8, 2)
seizure$x <- ifelse(seizure$visit == 0, 0, 1)
mm <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data=seizure, corstr="exch", family=poisson)
summary(mm)
mm <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data=seizure, subset=id!=49, corstr="exch", family=poisson)
summary(mm)
seiz2 <- subset(seizure, visit > 0)
seiz2$lbase <- log(rep(subset(seizure, visit == 0)$y, rep(4, 59)) / 4)
seiz2$lage <- log(seiz2$age)
seiz2$v4 <- ifelse(seiz2$visit == 4, 1, 0)
m2 <- geese(y ~ lbase + trt + lbase:trt + age + v4, id = id,
            data=seiz2, subset=id!=49, corstr="exch", family=poisson)
}
\keyword{nonlinear}
\keyword{models}
