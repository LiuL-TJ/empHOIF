## The following code is used to estimate the nuisance functions by lm and density g by kde
## Input: Daubechies_***.RData, the simulated dataset
## Output: Highdim_Nuisance_Ln_***.RData, including the nuisance function estimators by lm and hat(g) by kde

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
## subset_index <- as.numeric(args[13])

resolution <- 15

set.seed(1024)

## M.n <- round(log2(k.n - 4))
N <- 10000
grps <- ceiling(ntr / N)

load(paste0('Daubechies_', filter.number, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', 50, '.RData'))
Xtrain <- PStrain <- ORtrain <- Ytrain <- Atrain <- NULL
for (i in 1:d) {
    Xtrain <- cbind(Xtrain, as.numeric(sapply(1:grps, function (j) as.matrix(Sim.Data[[i + (j - 1) * d]][[1]])))[1:ntr])
    PStrain <- cbind(PStrain, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[5]]))[1:ntr])
    ORtrain <- cbind(ORtrain, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[6]]))[1:ntr])
}
ORcoeffs <- runif(d, -0.5, 0.5)
PScoeffs <- runif(d, -0.5, 0.5)
PS_train <- expit(PStrain %*% PScoeffs)
OR_train <- ORtrain %*% ORcoeffs
Y_train <- rnorm(ntr, mean = OR_train, sd = varY^.5)
A_train <- rbinom(ntr, size = 1, prob = PS_train)
Data_train <- as.data.frame(cbind(Xtrain, A_train, Y_train, PS_train, OR_train))
colnames(Data_train) <- c(sapply(1:d, function (j) paste0('X', j)), 'A', 'Y', 'PS', 'OR')

bandw <- Hscv(Data_train[, 1:d])
if (d <= 4) {
    ghat <- kde(x = Data_train[, 1:d], H = bandw, binned = TRUE)
} else {
    ghat <- kde(x = Data_train[, 1:d], H = bandw, binned = FALSE)
}
bandwA <- Hscv(Data_train[Data_train$A == 1, 1:d])
if (d <= 4) {
    ghatA <- kde(x = Data_train[Data_train$A == 1, 1:d], h = bandwA, binned = TRUE)
} else {
    ghatA <- kde(x = Data_train[Data_train$A == 1, 1:d], h = bandwA, binned = FALSE)
}

## bandw <- Hscv(Data_train[, 1:d])
## ghat <- kde(x = Data_train[, 1:d], H = bandw, binned = TRUE)

## save(bandw, ghat, file = paste0('MultiDim_Density_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', filter.number.analysis, '_', n, '_', d, '.RData'))

ORmod <- lm(formula = as.formula(paste0('Y ~ ', paste0(sapply(1:d, function (j) paste0('X', j)), collapse = " + "))), data = Data_train)
PSmod <- glm(formula = as.formula(paste0('A ~ ', paste0(sapply(1:d, function (j) paste0('X', j)), collapse = " + "))), data = Data_train, family = binomial)

save(Data_train, ORcoeffs, PScoeffs, ORmod, PSmod, bandw, ghat, bandwA, ghatA, file = paste0('Highdim_Nuisance_Ln_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))
