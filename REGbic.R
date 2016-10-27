#'  REGbic 
#'
#' This calulates regression BIC
#' @param y is dependent variable as numeric vector, and x is independent variable/s  a numeric matrix
#' @return  regression BIC (numeric)
#' @export
#' @examples
#' REGbic()


REGbic <- function (y, x) {

  y <- as.vector(y)
  x <- as.matrix(x)
  
  if (any(is.na(y))|| any(is.na(x))) {
    warning("NA's in the y or x")
    return(NULL)
  }
  
  if (any(is.null(y))|| any(is.null(x))) {
    return(NULL)
  }
  
  p <- ncol(x) + 2

  n  <- length(y)
  fit <- lm(y~x)
  sigma <- (sum((summary(fit)$resid)^2)/n)^0.5

  if (ncol(x)==1){
    REG.bic <- - (-n*log(2*pi)-2*n*log(sigma)-n-log(n)*3)
  } else {

    REG.bic <- - (-n*log(2*pi)-2*n*log(sigma)-n-log(n)*p)
  }
  return(REG.bic)
}

