rm(list = ls())
ptm <- proc.time()
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
library(Rfast)
library(mgcv)

args = commandArgs(trailingOnly = TRUE)

expit <- function (x) exp(x) / (1 + exp(x))

setwd('~/')
source('HOIF_Functions_Mat.R')

new.seed <- as.numeric(args[1])
beta.b <- as.numeric(args[2])
beta.p <- as.numeric(args[3])
beta.g <- as.numeric(args[4])
Rb <- as.numeric(args[5])
Rp <- as.numeric(args[6])
Rg <- as.numeric(args[7])
varY <- as.numeric(args[8])
filter.number <- as.numeric(args[9])
filter.number.analysis <- as.numeric(args[10])
ntr <- as.numeric(args[11])
d <- as.numeric(args[12])
n <- as.numeric(args[13])
k.n <- as.numeric(args[14])

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
load(paste0('Highdim_Nuisance_mGCV_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', ntr, '_', d, '.RData'))
load(paste0('Cov_Highdim_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', filter.number.analysis, '_', ntr, '_', d, '_', k.n, '.RData'))

load(paste0('Daubechies_', filter.number, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', new.seed, '.RData'))

grps <- ceiling(n / N)
X <- PS <- OR <- Y <- A <- NULL
for (i in 1:d) {
    X <- cbind(X, as.numeric(sapply(1:grps, function (j) as.matrix(Sim.Data[[i + (j - 1) * d]][[1]])))[1:n])
    PS <- cbind(PS, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[5]]))[1:n])
    OR <- cbind(OR, as.numeric(sapply(1:grps, function (j) Sim.Data[[i + (j - 1) * d]][[6]]))[1:n])
}
PS_true <- expit(PS %*% PScoeffs)
OR_true <- OR %*% ORcoeffs
Y <- rnorm(n, mean = OR_true, sd = varY^.5)
A <- rbinom(n, size = 1, prob = PS_true)
Data <- as.data.frame(cbind(X, A, Y, PS_true, OR_true))
colnames(Data) <- c(sapply(1:d, function (j) paste0('X', j)), 'A', 'Y', 'PS', 'OR')
Data$PS_est <- predict(PSgam, newdata = Data, type = 'response')
Data$OR_est <- predict(ORgam, newdata = Data)

Data$OR_res <- Data$Y - Data$OR_est
Data$PS_res <- 1 - Data$A / Data$PS_est

DB <- NULL
for (i in 1:d) {
    X <- as.numeric(Data_train[, i])
    DB_i <- wv_analysis(as.numeric(Data[, i]), wv.father, filter.number.analysis, M.n, 1:k.n)
    DB <- cbind(DB, DB_i)
}
DB <- DB[, filter.cols]

ptm2 <- proc.time()
Omega_train.sqrt <- mat.mult(mat.mult(svd.sigma.tr$v, diag(svd.sigma.tr$d^(-1/2))), t(svd.sigma.tr$u))
left.basis <- mat.mult(DB, Omega_train.sqrt)

psi.mar <- mean(Data$Y)
IF1.mar.mis <- with(Data, mean(OR_est + A / PS_est * OR_res))
IF1.mar.null <- with(Data, mean(0 + A / 0.5 * Y))
IF22.mar.emp <- HOIF_2_fast(Data$OR_res * Data$A, Data$PS_res, left.basis)
Data$PS_null <- 0.5
Data$PS_res_null <- Data$A - Data$PS_null
IF22.mar.emp.null <- HOIF_2_fast(Data$Y * Data$A, Data$PS_res_null, left.basis)

library(ks)
Data$fhat <- predict(ghatA, x = Data[, 1:d])
left.basis <- DB / (Data$fhat^(0.5) * mean(Data_train$A)^(0.5))
IF22.mar.ker <- HOIF_2_fast(Data$OR_res * Data$A, Data$PS_res, left.basis)
IF22.mar.ker.null <- HOIF_2_fast(Data$Y * Data$A, Data$PS_res_null, left.basis)

load(paste0('Highdim_Nuisance_AC_', filter.number, '_', filter.number.analysis, '_', varY, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', ntr, '_', d, '_', k.n, '.RData'))
OmegaA <- OmegaA * 1000 / rr
OmegaAAC <- OmegaA[filter.cols, filter.cols] * mean(Data_train$A)
svd.omegaAac <- svd(OmegaAAC, LINPACK = TRUE)
PrecisionAAC.sqrt <- svd.omegaAac$v %*% diag(svd.omegaAac$d^(-1/2)) %*% t(svd.omegaAac$u)
left.basis <- DB %*% PrecisionAAC.sqrt
IF22.mar.ac <- HOIF_2_fast(Data$OR_res * Data$A, Data$PS_res, left.basis)
IF22.mar.ac.null <- HOIF_2_fast(Data$Y * Data$A, Data$PS_res_null, left.basis)

load(paste0('Oracle_Higdim_', beta.p, '_', beta.g, '_', Rp, '_', Rg, '_', filter.number, '_', filter.number.analysis, '_', d, '_', k.n, '.RData'))
left.basis <- mat.mult(DB, Omega.sqrt)
IF22.mar.oracle <- HOIF_2_fast(Data$OR_res * Data$A, Data$PS_res, left.basis)
IF22.mar.oracle.null <- HOIF_2_fast(Data$Y * Data$A, Data$PS_res_null, left.basis)

save(psi.mar, IF1.mar.mis, IF1.mar.null, IF22.mar.emp, IF22.mar.emp.null, IF22.mar.ker, IF22.mar.ker.null, IF22.mar.ac, IF22.mar.ac.null, IF22.mar.oracle, IF22.mar.oracle.null, file = paste0('IF22_Highdim_mGCV_', new.seed, '_', beta.b, '_', beta.p, '_', beta.g, '_', Rb, '_', Rp, '_', Rg, '_', varY, '_', filter.number, '_', filter.number.analysis, '_', ntr, '_', d, '_', n, '_', k.n, '.RData'))
## save(psi.mar, IF1.mar.mis, IF1.mar.null, IF22.mar.emp, IF22.mar.emp.null, IF22.mar.ker, IF22.mar.ac, file = paste0('IF22_Highdim_Ln_', new.seed, '_', beta.b, '_', beta.p, '_', Rb, '_', Rp, '_', varY, '_', filter.number, '_', filter.number.analysis, '_', ntr, '_', d, '_', n, '_', k.n, '.RData'))
print(proc.time() - ptm)
print(proc.time() - ptm2)
