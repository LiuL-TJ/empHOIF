## Combining all the simulated true Omega together

rm(list = ls())

library(MASS)
library(wavethresh)
library(DaubechiesSim)
library(Rfast)

args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))

setwd('/net/rcstorenfs02/ifs/rc_labs/robins_lab/users/linliu824/HOIF/Empirical-HOIF/Revision/')

beta.p <- as.numeric(args[1])
beta.g <- as.numeric(args[2])
Rp <- as.numeric(args[3])
Rg <- as.numeric(args[4])
filter.number <- as.numeric(args[5])
filter.number.analysis <- as.numeric(args[6])
d <- as.numeric(args[7])
k.n <- as.numeric(args[8])

load(paste0('Cov_Highdim_', 0.25, '_', beta.p, '_', beta.g, '_', 1, '_', Rp, '_', Rg, '_', 1, '_', filter.number, '_', filter.number.analysis, '_', 2e+05, '_', d, '_', k.n, '.RData'))

Sigma_oracle <- 0
for (seed in 1:50) {
    load(paste0('Oracle_Higdim_', seed, '_', beta.p, '_', beta.g, '_', Rp, '_', Rg, '_', filter.number, '_', filter.number.analysis, '_', d, '_', k.n, '.RData'))
    Sigma_oracle <- Sigma_oracle + Sigma / 50
}

Sigma_oracle <- Sigma_oracle[filter.cols, filter.cols]
svd.sigma <- svd(Sigma_oracle, LINPACK = TRUE)
Omega.sqrt <- mat.mult(mat.mult(svd.sigma$v, diag(svd.sigma$d^(-1/2))), t(svd.sigma$u))

save(svd.sigma, Omega.sqrt, file = paste0('Oracle_Higdim_', beta.p, '_', beta.g, '_', Rp, '_', Rg, '_', filter.number, '_', filter.number.analysis, '_', d, '_', k.n, '.RData'))
