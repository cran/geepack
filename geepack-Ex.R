pkgname <- "geepack"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('geepack')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("DATA-dietox")
### * DATA-dietox

flush(stderr()); flush(stdout())

### Name: dietox
### Title: Growth curves of pigs in a 3x3 factorial experiment
### Aliases: dietox
### Keywords: datasets

### ** Examples

data(dietox)
dietox$Cu     <- as.factor(dietox$Cu)
gee01 <- geeglm (Weight ~ Time + Cu + Cu * Time, id =Pig, data = dietox,
         family=gaussian,corstr="ex")

mf <- formula(Weight~Cu*(Time+I(Time^2)+I(Time^3)))
gee1 <- geeglm(mf, data=dietox, id=Pig, family=poisson("identity"),corstr="ar1")
summary(gee1)
anova(gee1)




cleanEx()
nameEx("DATA-koch")
### * DATA-koch

flush(stderr()); flush(stdout())

### Name: koch
### Title: Ordinal Data from Koch
### Aliases: koch
### Keywords: datasets

### ** Examples

data(koch)
fit <- ordgee(ordered(y) ~ trt + as.factor(day), id=id, data=koch, corstr="exch")
summary(fit)



cleanEx()
nameEx("DATA-ohio")
### * DATA-ohio

flush(stderr()); flush(stdout())

### Name: ohio
### Title: Ohio Children Wheeze Status
### Aliases: ohio
### Keywords: datasets

### ** Examples

data(ohio)
fit <- geese(resp ~ age + smoke + age:smoke, id=id, data=ohio,
             family=binomial, corstr="exch", scale.fix=TRUE)
summary(fit)
fit.ar1 <- geese(resp ~ age + smoke + age:smoke, id=id, data=ohio,
                 family=binomial, corstr="ar1", scale.fix=TRUE)
summary(fit.ar1)



cleanEx()
nameEx("DATA-sitka89")
### * DATA-sitka89

flush(stderr()); flush(stdout())

### Name: sitka89
### Title: Growth of Sitka Spruce Trees
### Aliases: sitka89
### Keywords: datasets

### ** Examples

data(sitka89)



cleanEx()
nameEx("DATA-spruce")
### * DATA-spruce

flush(stderr()); flush(stdout())

### Name: spruce
### Title: Log-size of 79 Sitka spruce trees
### Aliases: spruce
### Keywords: datasets

### ** Examples

data(spruce)
spruce$contr <- ifelse(spruce$ozone=="enriched", 0, 1)
sitka88 <- spruce[spruce$wave <= 5,]
sitka89 <- spruce[spruce$wave > 5,]
fit.88 <- geese(logsize ~ as.factor(wave) + contr +
                          I(time/100*contr) - 1,
                id=id, data=sitka88, corstr="ar1")
summary(fit.88)

fit.89 <- geese(logsize ~ as.factor(wave) + contr - 1,
                id=id, data=sitka89, corstr="ar1")
summary(fit.89)



cleanEx()
nameEx("fixed2Zcor")
### * fixed2Zcor

flush(stderr()); flush(stdout())

### Name: fixed2Zcor
### Title: Construct zcor vector (of fixed correlations) from a fixed
###   working correlation matrix
### Aliases: fixed2Zcor
### Keywords: regression

### ** Examples


timeorder <- rep(1:5, 6)
tvar      <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each=5)
uuu   <- rep(rnorm(6), each=5)
yvar  <- 1 + 2*tvar + uuu + rnorm(length(tvar))
simdat <- data.frame(idvar, timeorder, tvar, yvar)
head(simdat,12)

library(doBy)
simdatPerm <- simdat[sample(nrow(simdat)),]
simdatPerm <- orderBy(~idvar, simdatPerm)
head(simdatPerm)

cor.fixed <- matrix(c(1    , 0.5  , 0.25,  0.125, 0.125,
                      0.5  , 1    , 0.25,  0.125, 0.125,
                      0.25 , 0.25 , 1   ,  0.5  , 0.125,
                      0.125, 0.125, 0.5  , 1    , 0.125,
                      0.125, 0.125, 0.125, 0.125, 1     ), 5, 5)
cor.fixed

zcor <- fixed2Zcor(cor.fixed, id=simdatPerm$idvar, waves=simdatPerm$timeorder)
zcor

mod4 <- geeglm(yvar~tvar, id=idvar, data=simdatPerm, corstr="fixed", zcor=zcor)
mod4





cleanEx()
nameEx("geeglm")
### * geeglm

flush(stderr()); flush(stdout())

### Name: geeglm
### Title: Fit Generalized Estimating Equations (GEE)
### Aliases: geeglm
### Keywords: models

### ** Examples

data(dietox)
dietox$Cu     <- as.factor(dietox$Cu)
mf <- formula(Weight~Cu*(Time+I(Time^2)+I(Time^3)))
gee1 <- geeglm(mf, data=dietox, id=Pig, family=poisson("identity"),corstr="ar1")
gee1
summary(gee1)

mf2 <- formula(Weight~Cu*Time+I(Time^2)+I(Time^3))
gee2 <- geeglm(mf2, data=dietox, id=Pig, family=poisson("identity"),corstr="ar1")
anova(gee2)




cleanEx()
nameEx("geese")
### * geese

flush(stderr()); flush(stdout())

### Name: geese
### Title: Function to solve a Generalized Estimating Equation Model
### Aliases: geese geese.fit print.geese summary.geese print.summary.geese
### Keywords: nonlinear models

### ** Examples

data(seizure)
## Diggle, Liang, and Zeger (1994) pp166-168, compare Table 8.10
seiz.l <- reshape(seizure,
                  varying=list(c("base","y1", "y2", "y3", "y4")),
                  v.names="y", times=0:4, direction="long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)
m1 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data=seiz.l, corstr="exch", family=poisson)
summary(m1)
m2 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data = seiz.l, subset = id!=49,
            corstr = "exch", family=poisson)
summary(m2)
## Using fixed correlation matrix
cor.fixed <- matrix(c(1, 0.5, 0.25, 0.125, 0.125,
                      0.5, 1, 0.25, 0.125, 0.125,
                      0.25, 0.25, 1, 0.5, 0.125,
                      0.125, 0.125, 0.5, 1, 0.125,
                      0.125, 0.125, 0.125, 0.125, 1), 5, 5)
cor.fixed
zcor <- rep(cor.fixed[lower.tri(cor.fixed)], 59)
m3 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data = seiz.l, family = poisson,
            corstr = "fixed", zcor = zcor)
summary(m3)

data(ohio)
fit <- geese(resp ~ age + smoke + age:smoke, id=id, data=ohio,
             family=binomial, corstr="exch", scale.fix=TRUE)
summary(fit)
fit.ar1 <- geese(resp ~ age + smoke + age:smoke, id=id, data=ohio,
                 family=binomial, corstr="ar1", scale.fix=TRUE)
summary(fit.ar1)

###### simulated data
## a function to generate a dataset
gendat <- function() {
  id <- gl(50, 4, 200)
  visit <- rep(1:4, 50)
  x1 <- rbinom(200, 1, 0.6) ## within cluster varying binary covariate
  x2 <- runif(200, 0, 1)   ## within cluster varying continuous covariate
  phi <- 1 + 2 * x1         ## true scale model
  ## the true correlation coefficient rho for an ar(1)
  ## correlation structure is 0.667.
  rhomat <- 0.667 ^ outer(1:4, 1:4, function(x, y) abs(x - y))
  chol.u <- chol(rhomat)
  noise <- as.vector(sapply(1:50, function(x) chol.u %*% rnorm(4)))
  e <- sqrt(phi) * noise
  y <- 1 + 3 * x1 - 2 * x2 + e
  dat <- data.frame(y, id, visit, x1, x2)
  dat
}

dat <- gendat()
fit <- geese(y ~ x1 + x2, id = id, data = dat, sformula = ~ x1,
             corstr = "ar1", jack = TRUE, j1s = TRUE, fij = TRUE)
summary(fit)


#### create user-defined design matrix of unstrctured correlation.
#### in this case, zcor has 4*3/2 = 6 columns, and 50 * 6 = 300 rows
zcor <- genZcor(clusz = rep(4, 50), waves = dat$visit, "unstr")
zfit <- geese(y ~ x1 + x2, id = id, data = dat, sformula = ~ x1,
              corstr = "userdefined", zcor = zcor,
              jack = TRUE, j1s = TRUE, fij = TRUE)
summary(zfit)

#### Now, suppose that we want the correlation of 1-2, 2-3, and 3-4
#### to be the same. Then zcor should have 4 columns.
z2 <- matrix(NA, 300, 4)
z2[,1] <- zcor[,1] + zcor[,4] + zcor[,6]
z2[,2:4] <- zcor[, c(2, 3, 5)]
summary(geese(y ~ x1 + x2, id = id, data = dat, sformula = ~ x1,
              corstr = "userdefined", zcor = z2,
              jack = TRUE, j1s = TRUE, fij = TRUE))

#### Next, we introduce non-constant cluster sizes by
#### randomly selecting 60 percent of the data
good <- sort(sample(1:nrow(dat), .6 * nrow(dat))) 
mdat <- dat[good,]

summary(geese(y ~ x1 + x2, id = id, data = mdat, waves = visit,
              sformula = ~ x1, corstr="ar1",
              jack = TRUE, j1s = TRUE, fij = TRUE))




cleanEx()
nameEx("genZcor")
### * genZcor

flush(stderr()); flush(stdout())

### Name: genZcor
### Title: genZcor
### Aliases: genZcor humbelbee
### Keywords: regression

### ** Examples

#example to construct a Toeplitz correlation structure
#    sigma_ij=sigma_|i-j|

#data set with 5 clusters and maximally 4 observations (visits) per cluster
 gendat <- function() {
       id <- gl(5, 4, 20)
       visit <- rep(1:4, 5)
       y <- rnorm(id)
       dat <- data.frame(y, id, visit)[c(-2,-9),]
}

set.seed(88)
dat<-gendat()

#generating the design matrix for the unstructured correlation
zcor <- genZcor(clusz = table(dat$id), waves = dat$visit, corstrv=4)
# defining the Toeplitz structure 
zcor.toep<-matrix(NA, nrow(zcor),3)
zcor.toep[,1]<-apply(zcor[,c(1,4,6)],1,sum)
zcor.toep[,2]<-apply(zcor[,c(2,5)],1,sum)
zcor.toep[,3]<-zcor[,3]

zfit1 <- geese(y ~ 1,id = id, data = dat,
                   corstr = "userdefined", zcor = zcor.toep)


zfit2 <- geeglm(y ~ 1,id = id, data = dat,
                   corstr = "userdefined", zcor = zcor.toep)





cleanEx()
nameEx("ordgee")
### * ordgee

flush(stderr()); flush(stdout())

### Name: ordgee
### Title: GEE for Clustered Ordinal Responses
### Aliases: ordgee
### Keywords: nonlinear models

### ** Examples

data(respdis)
resp.l <- reshape(respdis, varying =list(c("y1", "y2", "y3", "y4")),
                  v.names = "resp", direction = "long")
resp.l <- resp.l[order(resp.l$id, resp.l$time),]
fit <- ordgee(ordered(resp) ~ trt, id=id, data=resp.l, int.const=FALSE)
summary(fit)

data(ohio)
ohio$resp <- ordered(as.factor(ohio$resp))
fit <- ordgee(resp ~ age + smoke + age:smoke, id = id, data=ohio)
summary(fit)



cleanEx()
nameEx("respdis")
### * respdis

flush(stderr()); flush(stdout())

### Name: respdis
### Title: Clustered Ordinal Respiratory Disorder
### Aliases: respdis
### Keywords: datasets

### ** Examples

data(respdis)
resp.l <- reshape(respdis, varying = list(c("y1", "y2", "y3", "y4")),
                  v.names = "resp", direction = "long")
resp.l <- resp.l[order(resp.l$id, resp.l$time),]
fit <- ordgee(ordered(resp) ~ trt, id = id, data = resp.l, int.const = FALSE)
summary(fit)

z <- model.matrix( ~ trt - 1, data = respdis)
ind <- rep(1:111, 4*3/2 * 2^2)
zmat <- z[ind,,drop=FALSE]
fit <- ordgee(ordered(resp) ~ trt, id = id, data = resp.l, int.const = FALSE,
              z = zmat, corstr = "exchangeable")
summary(fit)



cleanEx()
nameEx("respiratory")
### * respiratory

flush(stderr()); flush(stdout())

### Name: respiratory
### Title: Data from a clinical trial comparing two treatments for a
###   respiratory illness
### Aliases: respiratory respiratoryWide
### Keywords: datasets

### ** Examples

data(respiratory)
## maybe str(respiratory) ; plot(respiratory) ...



cleanEx()
nameEx("seizure")
### * seizure

flush(stderr()); flush(stdout())

### Name: seizure
### Title: Epiliptic Seizures
### Aliases: seizure
### Keywords: datasets

### ** Examples

data(seizure)
## Diggle, Liang, and Zeger (1994) pp166-168, compare Table 8.10
seiz.l <- reshape(seizure,
                  varying=list(c("base","y1", "y2", "y3", "y4")),
                  v.names="y", times=0:4, direction="long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)
m1 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data=seiz.l, corstr="exch", family=poisson)
summary(m1)
m2 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
            data = seiz.l, subset = id!=49,
            corstr = "exch", family=poisson)
summary(m2)

## Thall and Vail (1990)
seiz.l <- reshape(seizure, varying=list(c("y1","y2","y3","y4")),
                  v.names="y", direction="long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time),]
seiz.l$lbase <- log(seiz.l$base / 4)
seiz.l$lage <- log(seiz.l$age)
seiz.l$v4 <- ifelse(seiz.l$time == 4, 1, 0)
m3 <- geese(y ~ lbase + trt + lbase:trt + lage + v4, 
            sformula = ~ as.factor(time) - 1, id = id,
            data = seiz.l, corstr = "exchangeable", family=poisson)
## compare to Model 13 in Table 4, noticeable difference
summary(m3)

## set up a design matrix for the correlation
z <- model.matrix(~ age, data = seizure)  # data is not seiz.l
## just to illustrate the scale link and correlation link
m4 <- geese(y ~ lbase + trt + lbase:trt + lage + v4,
            sformula = ~ as.factor(time)-1, id = id,
            data = seiz.l, corstr = "ar1", family = poisson,
            zcor = z, cor.link = "fisherz", sca.link = "log")
summary(m4)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
