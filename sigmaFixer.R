sigmaFixer <- function(sigma, ...) {
    sigma <- as.matrix(sigma)
    r <- dim(sigma)[1]
    fixrate <- 0.01
    covfixmat <- array(1, c(r, r)) + fixrate * diag(1, r)
    min_limit <- .Machine$double.eps * 10
    
    if (!all(is.finite(sigma))) {
        warning("covariance matrix is not finite")
    }
    
    # Enforces squareness and symmetricity
    nsigma <- sigma - ((sigma - t(sigma))/2)
    iter <- 0
    
    # Checking the covariance matrix is not positive definite
    while (postDef(nsigma) == 0 & (iter < 10000)) {
        
        iter <- iter + 1
        d <- diag(nsigma)
        if (any(d <= min_limit)) {
            m <- max(abs(d)) * fixrate
            neg <- min(d)
            
            if (neg < 0) {
                addit <- (m - neg) * diag(1, r)
            } else {
                if (m < min_limit) {
                  m <- min_limit
                }
                addit <- m * diag(1, r)
            }
            nsigma <- nsigma + addit
        } else {
            # Increase the diagonal values by 1%
            nsigma <- nsigma * covfixmat
        }
    }
    # return(list(nsigma, iter))
    return(nsigma)
}


# testing for postive definite
postDef <- function(Sigma) {
    Sigma <- as.matrix(Sigma)
    q <- try(chol(Sigma)[2], TRUE)
    if (q == 0) {
        t <- 1
    } else {
        t <- 0
    }
    return(t)
} 
