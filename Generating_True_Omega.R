## The following code is used to approximately evaluate the true Omega by generating L = 10^(7) samples
## Input: Daubechies_***.RData, including the datasets
## Output: Oracle_Higdim_*.RData, including Omega

rm(list = ls())

library(MASS)
library(wavethresh)
library(DaubechiesSim)
library(Rfast)

args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))

setwd('~/')

seed <- as.numeric(args[1])
beta.p <- as.numeric(args[2])
beta.g <- as.numeric(args[3])
Rp <- as.numeric(args[4])
Rg <- as.numeric(args[5])
filter.number <- as.numeric(args[6])
filter.number.analysis <- as.numeric(args[7])
d <- as.numeric(args[8])
k.n <- as.numeric(args[9])
resolution <- 15
M.n <- round(log2(k.n - 4))

wv.mother <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
wv.father <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)
wv.mother$x <- wv.mother$x + filter.number.analysis
haar.mother <- draw(filter.number = 1, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
haar.father <- draw(filter.number = 1, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)
haar.mother$x <- haar.mother$x + 1

res <- 6
X_space <- seq(0, 1, 10^(-res))

## load(paste0('Highdim_Nuisance_Ln_', 0.25, '_', beta.p, '_', beta.g, '_', 1, '_', Rp, '_', Rg, '_', 1, '_', filter.number, '_', 25000, '_', d, '.RData'))
load(paste0('Daubechies_', filter.number, '_', 1, '_', 0.25, '_', beta.p, '_', beta.g, '_', 1, '_', Rp, '_', Rg, '_', seed, '.RData'))
Xtrain <- PStrain <- ORtrain <- Ytrain <- Atrain <- NULL
grps <- ceiling(200000 / 10000)
for (i in 101:(100 + d)) {
    Xtrain <- cbind(Xtrain, as.numeric(sapply(1:grps, function (j) as.matrix(Sim.Data[[i + (j - 1) * d]][[1]]))))
    PStrain <- cbind(PStrain, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[5]])))
    ORtrain <- cbind(ORtrain, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[6]])))
}
PS_train <- expit(PStrain %*% PScoeffs)
A_train <- rbinom(length(PS_train), size = 1, prob = PS_train)

DB_train <- NULL
for (i in 1:d) {
    X <- as.numeric(Xtrain[, i])
    DB_train_i <- wv_analysis(X, wv.father, filter.number.analysis, M.n, 1:k.n)
    DB_train <- cbind(DB_train, DB_train_i)
}

Sigma_oracle <- 0
for (i in 1:500) {
    id <- (nrow(DB_train) / 500 * (i - 1) + 1):(nrow(DB_train) / 500 * i)
    Sigma_oracle <- Sigma_oracle + mat.mult(t(A_train[id] * DB_train[id, ]), A_train[id] * DB_train[id, ]) / nrow(DB_train)
    Sigma <- Sigma_oracle * 500 / i
    save(Sigma, file = paste0('Oracle_Higdim_', seed, '_', beta.p, '_', beta.g, '_', Rp, '_', Rg, '_', filter.number, '_', filter.number.analysis, '_', d, '_', k.n, '.RData'))
}
## save(Sigma_oracle, file = paste0('Oracle_Higdim_', seed, '_', beta.p, '_', beta.g, '_', Rp, '_', Rg, '_', filter.number, '_', filter.number.analysis, '_', d, '_', k.n, '.RData'))
