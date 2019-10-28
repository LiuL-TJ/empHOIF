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

beta.b <- as.numeric(args[1])
beta.p <- as.numeric(args[2])
beta.g <- as.numeric(args[3])
Rb <- as.numeric(args[4])
Rp <- as.numeric(args[5])
Rg <- as.numeric(args[6])
varY <- as.numeric(args[7])
filter.number <- as.numeric(args[8])
filter.number.analysis <- as.numeric(args[9])
ntr <- as.numeric(args[10])
d <- as.numeric(args[11])
k.n <- as.numeric(args[12])

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
