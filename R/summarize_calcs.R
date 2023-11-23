#' Summary of repeated measurements of signature associations 
#'
#' @description Over repeated measurements stored in lists of lists, the
#' function returns 0 if half of the values in each subsampling calculation is 0
#' for half of replication experiments. Otherwise, the function returns the average
#' from the experiments where at least half of the values are non-zero. 
#' 
#' @param all.estimates List of repeated correlation estimates between signatures. Required.
#'
#' @import dplyr
#' @import tidyverse
#' @import ggplot2
#'
#' @return returns a symmetric matrix of signature interactions, diag is 0. 
#' 
#' @export

summarize_calcs = function(all.estimates ) {

    simplify_by_exp = lapply(all.estimates, simplify2array)
    all_simplified = simplify2array(simplify_by_exp)
   
    ## learning some parameters

    sig.dims = dim( all_simplified ) [1]
    
    sig.lengths = sig.dims[1]
    
    active.squares = sig.lengths * (sig.lengths - 1) / 2

    mini.square.size = ceiling( sqrt( dim(all_simplified)[3] ) )
    smp.mat = all.estimates[[1]][[1]]

    exp.count.row = ceiling(sqrt(length(all.estimates)))

    exp.length = length(all.estimates)
    
    k = 0

    
    out.mat = matrix(0, ncol = sig.lengths, nrow = sig.lengths,
                     dimnames = list(rownames(smp.mat), colnames(smp.mat) ) )
    
    for (i in 1:(sig.lengths - 1)) {
        for (j in (i+1):sig.lengths) {
            ##            cat(i, j, "\n")

            exp.values = sapply(1:exp.length,  function(repi) {
                ## checks if half of the results are other than 0
                ij.sample.runs = all_simplified[i, j, , repi]
                if (sum(abs(ij.sample.runs) > 0 )  > (length(ij.sample.runs) / 2)) {
                    return(mean(ij.sample.runs))
                } else { return(NA)}

            } )

            if (sum(is.na(exp.values))  >= (exp.length / 2) ) {
                out.mat[i, j] = 0
                out.mat[j, i] = 0 
            } else {
                out.mat[i, j] = mean(exp.values, na.rm = TRUE)
                out.mat[j, i] = out.mat[i,j]
            } 
        }
    }
    return(out.mat)
}

