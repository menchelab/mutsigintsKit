#' Pairwise signature cooccurrences across samples
#'
#' @details For each pair of signatures cooccurrence is calculated using
#' Fisher's exact test. Before the test is performed, 1 is added to the table.
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0. default: TRUE
#' @param ... arguments passed to cor.test
#'
#' @param maxval The value to which Inf values will be converted to.
#'
#' @return Signature cooccurrence matrix. Values correspond to Fisher's test
#' odds ratios.
#'
#' @export


cooccurrence = function(x, p.val = 0.05, maxval = 10, p.adjust = TRUE,  ...) {

    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    else {
        x[ x > 0 ] = 1 
    }
    if (any(x[!is.na(x)] < 0)) 
        stop("all elements of x must be non-negative")

    cooc_matrix = matrix(0, nrow = ncol(x), ncol = ncol(x),
                         dimnames = list(colnames(x), colnames(x)))

    pvals_matrix = matrix(1, nrow = ncol(x), ncol = ncol(x),
                         dimnames = list(colnames(x), colnames(x)))

    for(i in 1:(ncol(cooc_matrix) - 1) ) {
        for(j in (i+1):ncol(cooc_matrix) ) {
            occurrences = table(x[, c(i, j) ] ) 
            if (sum(dim(occurrences)) < 4) next # pairs where one of the p(sig) == 1
            occ_test = fisher.test(occurrences)
            if (occ_test$p.value < p.val) {

                pvals_matrix[i,j] = occ_test$p.value # saving for the adjustement
                
                corr_status =  log(occ_test$estimate) ##
                if (is.infinite(corr_status)) {
                    corr_status = sign(corr_status) * maxval
                }
                ## corr_status =  occ_test$estimate

                cooc_matrix[i, j] = corr_status
            }
        }
    }

    if (p.adjust == TRUE) {

        calculated_ps = pvals_matrix[upper.tri(pvals_matrix)]
        adjusted_ps = p.adjust(calculated_ps, method = "BH")
        cooc_matrix[upper.tri(cooc_matrix)][which(adjusted_ps > p.val) ] = 0
    }


    cooc_matrix = cooc_matrix + t(cooc_matrix) ## mirror the upper.tri matrix
    
    omnipresent_sigs = names(which(colMeans(x) == 1 ) )


    if (length(omnipresent_sigs) == 1) {
        cooc_matrix[ omnipresent_sigs, omnipresent_sigs] = 1
    } else { 
        diag(cooc_matrix[ omnipresent_sigs, omnipresent_sigs]) = 1
    }

    return(cooc_matrix)
}
