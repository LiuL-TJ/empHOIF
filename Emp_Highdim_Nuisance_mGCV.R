## The following code is used to estimate the nuisance functions by GAM
## Input: Highdim_Nuisance_Ln_***.RData, including the training datasets
## Output: Highdim_Nuisance_mGCV_***.RData, including the nuisance function estimators by GAM

rm(list = ls())

library(MASS)
library(clusterGeneration)
library(wavethresh)
library(orthopolynom)
library(fda)
library(splines)
library(Rcpp)
library(parallel)
library(doMC)
library(iterpc)
library(gtools)
library(iterators)
library(corpcor)
library(DaubechiesSim)
library(np)
library(ks)
library(mgcv)
## use_condaenv("r-tensorflow")

args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))
logit <- function (x) log(x / (1 - x))

setwd('~/')

## new.seed <- as.numeric(args[1])
beta.b <- as.numeric(args[1])
beta.p <- as.numeric(args[2])
beta.g <- as.numeric(args[3])
Rb <- as.numeric(args[4])
Rp <- as.numeric(args[5])
Rg <- as.numeric(args[6])
varY <- as.numeric(args[7])
filter.number <- as.numeric(args[8])
ntr <- as.numeric(args[9])
d <- as.numeric(args[10])
## k.n <- as.numeric(args[12])
## M.n <- log2(k.n - 4)
## subset_index <- as.numeric(args[13])

resolution <- 15

set.seed(1024)

load(paste0('Highdim_Nuisance_Ln_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))

ORgam <- gam(formula = as.formula(paste0('Y ~ ', paste0(sapply(1:d, function (j) paste0('s(X', j, ')')), collapse = " + "))), data = Data_train, family = gaussian)
PSgam <- gam(formula = as.formula(paste0('A ~ ', paste0(sapply(1:d, function (j) paste0('s(X', j, ')')), collapse = " + "))), data = Data_train, family = binomial)

save(ORgam, PSgam, file = paste0('Highdim_Nuisance_mGCV_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))
