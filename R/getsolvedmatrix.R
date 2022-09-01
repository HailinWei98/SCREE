#' function definitions ##### Return liner regression results of perturbation score

getsolvedmatrix <- function(Xm, Ym, lambda = 0.01) {
    # Amat=solve(Xmat,Ymat) # solve AX=B, or Xmat * A =Ymat
    TMmat_g = crossprod(Xm, Xm) + lambda * diag(ncol(Xm))

    Amat_g = Ym %*% Xm %*% solve(TMmat_g)
    return(Amat_g)
}
