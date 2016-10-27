#'  EMvSel 
#'
#' This function selects clustering relevant variables for model-based clustering using greedy forward selection algorithm 
#' @param X is a numeric matrix/dataframe, bic is a criteria to select a variable as clustering relevant, default is 0
#' @return  clustering relevant variables (character)
#'  @export
#' @examples
#' EMvSel()
########################################################################################## 


# forward feature selection algorithm Forward greedy search
EMvSel <- function(X, bic = 0) {
    
    X <- as.data.frame(X)
    
    if (is.null(X)) 
        return(NULL)
    
    # X is non-missi
    if (any(is.na(X))) {
        warning("NA's in the dataframe")
        return(NULL)
    }
    
    n <- nrow(X)
    d <- ncol(X)
    
    if (is.null(X)) 
        return(NULL)
    
    # First Step - selecting single variable
    
    BIC_opt <- rep(NA, d)  # optimal BIC
    BIC_diff <- rep(NA, d)
    BIC_one <- rep(NA, d)
    
    
   
    
    
    # identify the optimal univariate clusters
    
    univ.BICs <- lapply(2:3, function(i) try(apply(X, 2, function(X) gmmEM(X, 
        c = i, initialize = c("kmeans"))$BIC), TRUE))
    
    univ.BIC <- data.frame(t(sapply(univ.BICs, c)))
    
    univ.BIC <- replace(univ.BIC, is.na(univ.BIC), .Machine$double.xmax)
    
    try(BIC_one <- apply(X, 2, function(X) gmmEM(X, c = 1, initialize = c("kmeans"))$BIC), 
        TRUE)
    
    BIC_sd <- sweep(-(univ.BIC), 2, FUN = "+", BIC_one)  #BIC_d <- c(BIC_one - BIC_opt)
    
    BIC_d <- apply(BIC_sd, 2, max)  # choosing the highest difference 
    
    # Find the variable with the highest BIC difference between optimal
    # clusters and no cluster
    
    v <- max(BIC_d[is.finite(BIC_d)])
    g <- which(BIC_d == v, arr.ind = TRUE)[1]
    
    # This is the first selected variable with most univariate clustering
    # evidence
    if (max(BIC_d[is.finite(BIC_d)]) > 0) {
        S <- matrix(c(X[, g]), n, 1)
        
    } else {
        return(list(VarSel = NULL, var = "No relevant variable"))
        stop("No relevant variable")
    }
    
    # optimal BIC of the first selected relevant variable
    BIC_S <- min(univ.BIC[g])
    colnames(S) <- colnames(X)[g]
    
    
    
    if (ncol(X) == 1) {
        
        if (BIC_d[g] >= 10) {
            
            return(list(VarSel = colnames(X)))
        }
        if (BIC_d[g] < 10) {
            # return(NULL)
            return(list(VarSel = NULL, var = "No relevant variable"))
            stop("No relevant variable")
        }
    }
    
    
    
    if (ncol(X) > 1) {
        
        # N is the matrix of currently irrelevant variables
        N <- as.matrix(X[, -g])
        colnames(N) <- colnames(X)[-g]
        
        # subset is a matrix records the proposed variable, optimal BIC of S and
        # difference in BIC for clustering versus no clustering on S.
        
        subset <- matrix(c(colnames(S), round(BIC_S, 4), round(BIC_d[g], 
            4), "Yes"), 1, 4)
        
        
        # Second Step - selecting second variable
        BIC_reg <- rep(0, ncol(N))
        BIC_joint <- rep(0, ncol(N))
        BIC_sum <- rep(0, ncol(N))
        BIC_j <- rep(0, ncol(N))
        
        
        # Bivariate joint clustering
        
        biv.BICs <- lapply(2:6, function(i) try(apply(N, 2, function(N) gmmEM(cbind(S, 
            N), c = i, initialize = c("kmeans"))$BIC), TRUE))
        err <- sapply(biv.BICs, is, class2 = "try-error")
        nulls <- sapply(biv.BICs, is, class2 = "NULL")
        biv.BIC_s <- biv.BICs[err == FALSE & nulls == FALSE]
        biv.BIC <- data.frame(t(sapply(biv.BIC_s, c)))
        
        biv.BIC <- replace(biv.BIC, is.na(biv.BIC), .Machine$double.xmax)
        
        
        # regressing non-relevant variable on relevant variable
        try(BIC_reg <- apply(N, 2, function(N) REGbic(N, S)), TRUE)
        
        
        BIC_sum <- BIC_reg + BIC_S
        
        BIC.df <- sweep(-(biv.BIC), 2, FUN = "+", BIC_sum)
        
        
        BIC_diff <- apply(BIC.df, 2, max)
        
        
        
        # Choose the variable with the largest difference
        v <- max(BIC_diff[is.finite(BIC_diff)])
        g <- which(BIC_diff == v, arr.ind = TRUE)[1]
        
        BIC_opt <- apply(biv.BIC, 2, min)
        
        # add the second best variable if its BIC difference is positive
        if (BIC_diff[g] > bic) {
            subset <- rbind(subset, c(colnames(N)[g], round(BIC_opt[g], 
                4), round(BIC_diff[g], 4), "Yes"))
            j <- c(colnames(S), colnames(N)[g])
            S <- cbind(S, N[, g])
            colnames(S) <- j
            N <- as.matrix(N[, -g])
        } else {
            
            subset <- rbind(subset, c(colnames(N)[g], round(BIC_opt[g], 
                4), round(BIC_diff[g], 4), "No"))
            
            return(list(VarSel = colnames(S), Steps = subset))
            stop()
        }
        
        
        
        S <- as.data.frame(S)
        ss <- names(S)
        N <- X[, -which(names(X) %in% ss)]
        N <- as.data.frame(N)
        colnames(N) <- names(X)[!names(X) %in% names(S)]
        iter <- 0
        
        
        while ((ncol(N) != 0) & !is.null(ncol(N)) & (iter < ncol(X))) {
            
            iter <- iter + 1
            
            BIC_reg <- rep(0, ncol(N))
            BIC_joint <- rep(0, ncol(N))
            BIC_diff <- rep(0, ncol(N))
            BIC_opt <- rep(0, ncol(N))
            
            
            # identifying optimal multivariate clusters
            
            mv.BICs <- lapply(2:6, function(i) try(gmmEM(S[, j], c = i, 
                initialize = c("kmeans"))$BIC, TRUE))  #  optimal BIC using relevant variables
            
            mv.BIC <- data.frame(t(sapply(mv.BICs, c)))
            mv.BIC <- replace(mv.BIC, is.na(mv.BIC), .Machine$double.xmax)
            
            mr <- min(sapply(mv.BIC, cbind))  #   BIC for optimal cluster
            MC <- which(sapply(mv.BIC, cbind) == mr, arr.ind = TRUE)[1]
            MC <- MC + 1  #  optimal clusters 
            
            
            # BIC of optimal cluster
            try(BIC_opt <- gmmEM(S[, j], c = MC, initialize = c("kmeans"))$BIC, 
                TRUE)
            try(BIC_joint <- apply(N, 2, function(N) gmmEM(cbind(S[, j], 
                N), c = MC, initialize = c("kmeans"))$BIC), TRUE)
            
            
            # regressing non-relevant variable on relevant variable
            try(BIC_reg <- apply(N, 2, function(N) REGbic(N, S[, j])), TRUE)
            
            
            BIC_sum <- BIC_reg + BIC_opt
            BIC_diff <- BIC_sum - BIC_joint
            
            
            # Choose the variable with the largest difference
            v <- max(BIC_diff[is.finite(BIC_diff)])
            g <- which(BIC_diff == v, arr.ind = TRUE)[1]
            
            
            if (BIC_diff[g] > bic) {
                
                # if this difference is positive add this variable to S and update the
                # clustering model's BICs
                
                subset <- rbind(subset, c(colnames(N)[g], round(BIC_opt, 
                  4), round(BIC_diff[g], 4), "Yes"))
                j <- c(colnames(S), colnames(N)[g])
                S <- as.data.frame(cbind(S, N[, g]))
                colnames(S) <- j
                
                ss <- names(S)
                
                N <- (X[, -which(names(X) %in% ss)])
                N <- as.data.frame(N)
                colnames(N) <- names(X)[!names(X) %in% names(S)]
                
            } else {
                
                
                subset <- rbind(subset, c(colnames(N)[g], round(BIC_opt, 
                  4), round(BIC_diff[g], 4), "No"))
                break
            }
        }
        
        
        # List the selected variables and the matrix of steps' information
        colnames(subset) <- c("ProposedVariable", "optimal_BIC", "BIC_Difference", 
            "Relevant?")
        return(list(VarSel = colnames(S), Steps = subset))
    }
} 
