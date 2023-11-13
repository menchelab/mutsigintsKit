#' Pairwise signature correlations from symmetric pivot coordinates
#'
#' @details For each pair of signatures correlation based on symmetric pivot
#' coordinates. The non-CoDa version of the function is \code{\link{cor_sigs}}.
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param rand.add If rand.add = TRUE, then the matrix is strictly positive,
#' otherwise rand.add = FALSE(default). Check what you're  
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0.
#' @param ... arguments passed to cor.test
#'
#' @return Signature correlation matrix
#'
#' @import mpmi
#' 
#' @export


cor_coda = function(x,  p.val = 0.05, rand.add = FALSE,
                    p.adjust = TRUE, mi = FALSE,  ...) {

    if (any(x[!is.na(x)] <= 0)) {
        if (rand.add){
            stop("all elements of x must be greater than 0")
        }
    }
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
    if (ncol(x) <= 2) 
        stop("Calculation of average symmetric coordinates not possible. 2 or less columns provided.")

    ind <- c(1:ncol(x))
    corZav <- matrix(NA, ncol(x), ncol(x))
    corPvals <- matrix(NA, ncol(x), ncol(x))
    
    for (i in 1:(ncol(x) - 1)) {
        for (j in (i + 1):ncol(x)) {

            ## cat (i, j, "\n")
            
            permuted_x =  x[, c(i, j, ind[-c(i, j)])]

            
            if (nrow(permuted_x) < 3 ) {
                corPvals[i, j]  <-  1
                corZav[i, j] <- 0
                next
            }


            if (rand.add) {
                
                balZavout = balZav(permuted_x)

                ## balZavout = balZavout[is.finite(rowSums(balZavout)), ]
                
            } else {
            
                ## If random noise was not added, then symm_coords
                ## are calculated for samples individually
                ## ATM if one of the pair is 0, both are returned as 0
                
                balZavout = apply( as.data.frame(permuted_x), MARGIN = 1,
                                  function(rowvec) {
                                      if (rowvec[1] == 0 | rowvec[2] == 0 ) return(c(0,0))
                                      rowvec = rowvec[ which (rowvec > 0) ]
                                      if (length(rowvec) == 2) return(c(0,0))
                                      return(balZavRow(rowvec))
                                  } ) %>% t()
                

                balZavout = balZavout[rowSums(abs(balZavout)) > 0, , drop = FALSE]

                balZavout = balZavout[is.finite(rowSums(balZavout)), , drop = FALSE]
                
                if (nrow(balZavout) < 3 ) {
                    corPvals[i, j]  <-  1
                    corZav[i, j] <- 0
                    next
                }
            }


            ## if (nrow(balZavout) < 3) {} 
            
            if (mi == TRUE) { ### check if the mpmi MI calculation should be used

                mi.out = mpmi::cmi(balZavout)

                corZav[i, j] <- mi.out$bcmi[1,2]
                corPvals[i, j] <- 2 * ( -abs( min(abs(mi.out$zvalues[1,2]), 100) ))
                
            } else { ### cor.test is used instead
                
                ## print(balZavout)
                corout = cor.test(balZavout[,1], balZavout[, 2], ...)
                
                corPvals[i, j]  <-  corout$p.value
                corZav[i, j] <- corout$estimate
            }
        }
    }

    if(p.adjust) {
        corPvals = p.adjust(corPvals, method = "BH")
    }
    
    corZav[corPvals > p.val ] = 0
    corZav[lower.tri(corZav)] <- t(corZav)[lower.tri(corZav)]
    diag(corZav) <- 1

    colnames(corZav) = colnames(x)
    rownames(corZav) = colnames(x)
    
    return(corZav)
}
