SoS <- function (v) {
    return(sum(v^2))
}

HOIF_2_fast <- function (Y.resid, A.resid, left.basis, right.basis) {
    n <- length(Y.resid)
    term1 <- as.numeric(mat.mult(t(Y.resid), left.basis))
    term2 <- as.numeric(mat.mult(t(A.resid), right.basis))
    term12 <- sum(term1 * term2) / (n * (n - 1))
    term3 <- apply(left.basis, 1, SoS)
    term4 <- sum(as.numeric(Y.resid) * term3 * as.numeric(A.resid)) / (n * (n - 1))
    return(- (term12 - term4))
}
