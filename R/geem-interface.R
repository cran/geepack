## Coerce between geeglm and geem objects
##

## as.geem  <- function(object) UseMethod("as.geem")

## as.geem.geeglm  <- function(object){
##     object2 <- object$call
##     object2[[1]] <- as.name("geem")
##     eval(object2)
## }

## as.geeglm <- function(object) UseMethod("as.geeglm")

## as.geem.geem  <- function(object){
##     object2 <- object$call
##     object2[[1]] <- as.name("geeglm")
##     eval(object2)
## }
