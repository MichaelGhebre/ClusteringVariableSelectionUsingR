#' gmmVarSel 
#'
#' This function implements the function (subEMvSel) in each cluster to search further relevant varibles .
#' @param Z is a numeric  matrix/ dataframe; BIC_diff is the criteria to select relevant variable in each cluster, default is 10
#' @return final clustering relevant  variables (character)
#' @export
#' @examples 
#' gmmVarSel(iris[,1:4])
#' [1] 'Petal.Length' 'Sepal.Width'  'Petal.Width' 


gmmVarSel <- function(Z, BIC_diff = 10) {
    
    Z <- as.data.frame(Z)
    
    s <- suppressWarnings(subEMvSel(Z))  #  relevant variables selected using algorithm 1
    
    if (ncol(as.matrix(Z[, !colnames(Z) %in% s])) >= 1 & length(s) > 0) {
        
        opt_clusterBIC <- lapply(2:3, function(i) try(gmmEM(Z[, s], c = i, 
            initialize = c("kmeans"))$BIC, TRUE))  #  optimal BIC using relevant variables
        opt_clusterBIC <- sapply(opt_clusterBIC, cbind)
        opt_clusterBIC <- replace(opt_clusterBIC, is.na(opt_clusterBIC), 
            .Machine$double.xmax)
        vv <- min(sapply(opt_clusterBIC, cbind))  #   BIC for optimal cluster
        G <- which(sapply(opt_clusterBIC, cbind) == vv, arr.ind = TRUE)[1]
        G <- G + 1  #  optimal clusters using relevant variables
        
        class <- gmmEM(Z[, s], c = G, initialize = c("kmeans"))$class  # class assigment from algorithm 1
        Z.new <- cbind(Z, class)
        c_var <- suppressWarnings(lapply(1:G, function(i) try(subEMvSel(Z[Z.new$class == 
            i, !colnames(Z) %in% s], BIC = BIC_diff), TRUE)))  #  select relevant in each cluster
        
        c_s <- as.vector(unique(do.call(rbind, as.list(sapply(c_var, cbind)))))  #  relevant variable in each cluster
        cs <- colnames(Z[, colnames(Z) %in% c_s, drop = F])
        all_s <- unique(as.vector(c(s, cs)))  #  all selected relvant variables from algorithm 1 and 2
        
        return(all_s)
    } else {
        return(s)
        
    }
    
} 
