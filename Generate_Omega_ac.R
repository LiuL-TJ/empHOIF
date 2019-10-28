## The following code is used to simulate L = 10000*1000 = 10^(7) samples from the estimated density g, in order to evaluate 
## the numerical integration hat(Omega)^(ac)
## Input: Highdim_Nuisance_Ln_***.RData, which includes the nuisance function estimators by linear model and the kernel density
##        estimator hat(g) using kde with bandwidth selected by cross-validation
## Output: Highdim_Nuisance_AC_***.RData, which includes hat(Omega)^(ac) for either the functional E[Y, A | X] or E[A p(X) Y]

rm(list = ls())

library(MASS)
library(wavethresh)
library(splines)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(doMC)
library(foreach)
library(doSNOW)
library(doParallel)
library(fda)
library(inline)
library(RInside)
library(RcppParallel)
library(iterpc)
library(gtools)
library(iterators)
library(corpcor)
library(orthopolynom)
library(DaubechiesSim)
library(ks)
args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))

setwd('~/')

beta.b <- as.numeric(args[1])
beta.p <- as.numeric(args[2])
beta.g <- as.numeric(args[3])
Rb <- as.numeric(args[4])
Rp <- as.numeric(args[5])
Rg <- as.numeric(args[6])
varY <- as.numeric(args[7])

N <- 10000

filter.number <- as.numeric(args[8])
filter.number.analysis <- as.numeric(args[9])
ntr <- as.numeric(args[10])
d <- as.numeric(args[11])
k.n <- as.numeric(args[12])
M.n <- log2(k.n - 4)

resolution <- 15
wv.mother <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
wv.father <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)
wv.mother$x <- wv.mother$x + filter.number.analysis

## gammas <- c(0, 3, 6, 9, 10, 16)

## wv.mother <- draw(filter.number = filter.number, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
## wv.father <- draw(filter.number = filter.number, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)

## number_cores <- min(10, detectCores())
## cl <- makeCluster(number_cores - 1)
## registerDoParallel(cl)

Sim.Data <- list()

set.seed(1)

load(paste0('Highdim_Nuisance_Ln_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))

## X_space_g_density <- predict(ghat, x = X_spaces)
## XA_space_g_density <- predict(ghatA, x = X_spaces)
rep <- 1000
ptm <- proc.time()
Omega <- OmegaA <- 0

if (file.exists(paste0('Highdim_Nuisance_AC_', filter.number, '_', filter.number.analysis, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', ntr, '_', d, '_', k.n, '.RData'))) {
    load(paste0('Highdim_Nuisance_AC_', filter.number, '_', filter.number.analysis, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', ntr, '_', d, '_', k.n, '.RData'))
    rrr <- rr
} else {
    rrr <- 1
}	
for (rr in (rrr + 1):rep) {
    ## First step: generate the covariates and data
    ## X.id <- sample(1:nrow(X_spaces), size = N, prob = X_space_g_density)
    ## XA.id <- sample(1:nrow(X_spaces), size = N, prob = XA_space_g_density)
    ## X <- X_spaces[X.id, ]; XA <- X_spaces[XA.id, ]
    ## X_hist <- hist(X, breaks = 200, plot = FALSE)
    X <- rkde(fhat = ghat, n = 10000)
    XA <- rkde(fhat = ghatA, n = 10000)
    ## X <- pmax(pmin(X, 1-1e-3), 1e-3)
    ## XA <- pmax(pmin(XA, 1-1e-3), 1e-3)
    DB.cov <- DBA.cov <- NULL
    for (dd in 1:d) {
        Xi <- X[, dd]; XAi <- XA[, dd]
        Xi <- pmax(pmin(Xi, 1-1e-3), 1e-3)
        XAi <- pmax(pmin(XAi, 1-1e-3), 1e-3)
        DB.cov = cbind(DB.cov, wv_analysis(Xi, wv.father, filter.number.analysis, M.n, 1:k.n))
        DBA.cov = cbind(DBA.cov, wv_analysis(XAi, wv.father, filter.number.analysis, M.n, 1:k.n))
    }
    ## Second step: generate R and Y
    ## nuis_functions <- foreach(k = 1:length(X), .combine = rbind, .packages = c('DaubechiesSim')) %dopar% {
    ##    c(expit(wv_trans_fast(X[k], wv.mother, filter.number, Rp, gammas, beta.p)))
    ## }

    ## oregs <- foreach(k = 1:length(X), .combine = c, .packages = c('DaubechiesSim')) %dopar% {
    ##     wv_trans_fast(X[k], wv.mother, filter.number, Rb, gammas, beta.b)
    ## }
    ## pscores <- nuis_functions[, 1]

    ## R <- rbinom(N, 1, pmax(pscores, 0.1))
    
    Omega <- Omega + t(DB.cov) %*% DB.cov / (N * rep)
    OmegaA <- OmegaA + t(DBA.cov) %*% DBA.cov / (N * rep)
    save(Omega, OmegaA, rr, file = paste0('Highdim_Nuisance_AC_', filter.number, '_', filter.number.analysis, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', ntr, '_', d, '_', k.n, '.RData'))
}
## stopCluster(cl)
print(proc.time() - ptm)
