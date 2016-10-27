#' subEMvSel 
#'
#' This function splits the entire variables into independent small subsets, and implement the  function 'EMvSel' to each subset to select clustering  relevant variables
#' @param Y is a numeric  matrix/dataframe; BIC  is a criteria to include as relevant variable, default is 0
#' @return  aggregated clustering relevant variables from each variable subset (character)
#' @export
#' @examples
#' subEMvSel()
#' 
############################################################################################# blocks (subset) of variables were idenified using factor analysis with
############################################################################################# varimax rotation# The above function (EMvSel) algrithm is applid in
############################################################################################# each block and the relevant variables #


# Factor Analysis


subEMvSel <- function(Y, BIC = 0) {
    
    require(psych)
    Y <- as.data.frame(Y)
    d <- ncol(Y)
    
    # indentify the possible number of factors for the entire variables
    # using likelihood approach
    fact <- lapply(1:d, function(i) try(fa(Y, i, SMC = F, rotate = "varimax")$factors, 
        TRUE))
    
    f <- max(suppressWarnings(na.omit(as.numeric(fact))))
    
    if (f >= 2) {
        # if the variables have more than two blocks implement the following
        # codes
        loading <- abs(fa(Y, f, SMC = F, rotate = "varimax")$loadings[, 
            1:f])
        factors <- as.data.frame(loading)
        
        eigenvalue <- apply(factors, 2, function(x) sum(x^2))
        kk <- as.data.frame(eigenvalue)
        
        # extracting factor which have total explained variance greater than 1
        eigen1 <- kk[kk > 1, 1]
        
        k <- length(eigen1)  # possible number of factors
        
        if (k > 1) {
            loading <- fa(Y, k, SMC = F, rotate = "varimax")$loadings[, 
                1:k]
            loading <- abs(as.data.frame(loading))
            
            # identify subset of variables (creating block for variables before
            # clustering)
            subst <- apply(loading, 1, function(x) sample(c(colnames(loading)[which(x == 
                max(x))]), ))
            s1 <- as.vector(subst)
            s2 <- cbind(s1, variables = names(subst))
            # dt <- data.table(s2)
            dt <- as.data.frame(s2)
            dtt <- cbind(dt, Factor = as.numeric(dt$s1))
            dtt <- as.data.frame(dtt)
            datt <- dtt[, 2:3]
            
            f.subset <- as.character()
            q <- length(table(datt$Factor))  # final variables blocks
            
            for (i in 1:q) {
                sub <- datt[datt$Factor == i, ]
                rownames(sub) <- sub[, 1]
                subset <- rownames(sub)[rownames(sub) %in% names(Y)]
                mm <- EMvSel(Y[, c(subset), drop = F], bic = BIC)  # feature selection in each block and collect the relevant variables
                # s <- names(mm$VarSel)
                f.subset <- c(f.subset, c(mm$VarSel))
            }
            
            
        } else {
            f.subset <- EMvSel(Y, bic = BIC)$VarSel
            
        }
    } else {
        f.subset <- EMvSel(Y, bic = BIC)$VarSel  # extract these matching with the orginal variables only ...to avoid null
    }
    
    return(f.subset)
} 
