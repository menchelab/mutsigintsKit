#' Calculating bias-corrected mutual information for signature pairs
#'
#' @param x Data.frame or a matrix of signature intensities(or counts).
#' Signatures need to correspond to columns.
#'
#' @param p.val All the correlations with >=p.val are set to 0.
#' @param p.adjust If TRUE, p.value adjustment with BH correction is applied
#' before setting correlations with p.value >= p.val to 0.
#'
#' @return BCMI matrix
#'
#' @import mpmi
#' 
#' @export

bcmi = function(x,  p.val = 0.05, p.adjust = TRUE) {

    if (any(x[!is.na(x)] < 0)) {
        stop("all elements of x must be greater than or equal to 0")
    }
    
    if (!is.matrix(x) & !is.data.frame(x)) 
        stop("x must be a matrix or data.frame")
    
    if (ncol(x) < 2) 
        stop("Less than two columns provided. Execution stops.")

    ## cmi doesn't work for columns which have only 0-valued elements

    mci_out_mat = matrix(0, ncol = ncol(x), nrow = ncol(x))

    zero_indeces = which(colSums(x) == 0)
    calced_indeces = setdiff(1:ncol(x), zero_indeces)

    x_cl = x[, colSums(x) > 0]
    
    cmi_out = mpmi::cmi(x_cl)

    bcmi_out = cmi_out$bcmi
    zvals_out = cmi_out$zvalues

    pvals = pnorm(-abs(zvals_out) )

    if ( p.adjust) {
        pvals = p.adjust(pvals, method = "BH")
    }

    pvals = matrix(pvals, ncol = ncol(zvals_out))

    filtered_bcmi = bcmi_out
    filtered_bcmi[pvals > p.val] = 0

    mci_out_mat[calced_indeces, calced_indeces] = filtered_bcmi
    
    colnames(mci_out_mat) = colnames(x)
    rownames(mci_out_mat) = colnames(x)
    
    return(mci_out_mat)
}
