## The following code is used to evaluate the sample Gram matrix hat(Omega) from the training sample

rm(list = ls())

library(MASS)
library(wavethresh)
library(splines)
library(Rcpp)
library(RcppArmadillo)
library(fda)
library(RInside)
library(RcppParallel)
library(iterpc)
library(gtools)
library(iterators)
library(corpcor)
library(orthopolynom)
library(DaubechiesSim)

args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))

setwd('~/')
source('HOIF_Functions_Mat.R')

beta.b <- as.numeric(args[1]) ## Smoothness of outcome regression
beta.p <- as.numeric(args[2]) ## Smoothness of propensity score
beta.g <- as.numeric(args[3]) ## Smoothness of density of X
Rb <- as.numeric(args[4]) ## Constant for outcome regression: Rb = 1 in the simulation shown in the paper
Rp <- as.numeric(args[5]) ## Constant for outcome regression: Rp = -2 in the simulation shown in the paper
Rg <- as.numeric(args[6]) ## Constant for outcome regression: Rg = 0.5 in the simulation shown in the paper
varY <- as.numeric(args[7]) ## The homoscedastic conditional variance of Y|X
filter.number <- as.numeric(args[8]) ## The number of vanishing moments for the Daubechies wavelets used to generate the data 
filter.number.analysis <- as.numeric(args[9]) ## The number of vanishing moments for the Daubechies wavelets used to analyze the data
ntr <- as.numeric(args[10]) ## Number of training samples
d <- as.numeric(args[11]) ## Dimension of the covariates X
k.n <- as.numeric(args[12]) ## The number of basis used to construct the second order influence function

M.n <- log2(k.n - 4)

resolution <- 15

set.seed(1024)
N <- 10000

wv.mother <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
wv.father <- draw(filter.number = filter.number.analysis, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)
wv.mother$x <- wv.mother$x + filter.number.analysis
haar.mother <- draw(filter.number = 1, family = "DaubExPhase", resolution = 2^resolution, scaling.function = FALSE, plot.it = FALSE, enhance = FALSE)
haar.father <- draw(filter.number = 1, family = "DaubExPhase", resolution = 2^resolution, scaling.function = TRUE, plot.it = FALSE, enhance = FALSE)
haar.mother$x <- haar.mother$x + 1

full.ests <- NULL

j <- 1

load(paste0('Highdim_Nuisance_Ln_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))

DB_train <- NULL
filter.cols <- NULL
for (i in 1:d) {
    X <- as.numeric(Data_train[, i])
    DB_train_i <- wv_analysis(X, wv.father, filter.number.analysis, M.n, 1:k.n)
    filter.col <- which(colMeans(DB_train_i) != 0)
    if (is.null(DB_train) == TRUE) {
        filter.cols <- c(filter.cols, filter.col)
    } else {
        filter.cols <- c(filter.cols, filter.col + ncol(DB_train))
    }
    DB_train <- cbind(DB_train, DB_train_i[, filter.col])
}
## Sigma_train <- t(Data_train$A * DB_train) %*% (Data_train$A * DB_train) / nrow(DB_train)
## Sigma_train <- crossprod(t(Data_train$A * DB_train)) / nrow(DB_train)
Sigma_train <- 0
for (i in 1:100) {
    id <- (nrow(DB_train) / 100 * (i - 1) + 1):(nrow(DB_train) / 100 * i)
    Sigma_train <- Sigma_train + mat.mult(t(Data_train$A[id] * DB_train[id, ]), Data_train$A[id] * DB_train[id, ]) / nrow(DB_train)
}
svd.sigma.tr <- svd(Sigma_train, LINPACK = TRUE)

save(svd.sigma.tr, filter.cols, file = paste0('Cov_Highdim_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', filter.number.analysis, '_', ntr, '_', d, '_', k.n, '.RData'))
