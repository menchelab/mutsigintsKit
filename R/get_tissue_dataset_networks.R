#' Extracts the interaction metrics of a tissue
#' from a list of interaction networks.
#' @param tissue Tissue to be extracted
#' @param network.lists The input list
#' @param filter.list A list of metric values to be filtered out. The list
#' element names have the same names as metric names in network.lists.
#' The values specify filtering threshold. Everything below the number is
#' set to 0. Default: NULL
#' @param filter.mat If provided, this is a list of matrices, where the name of
#' the element is the metric name. The list elements are matrices of all
#' signature combinations in that tissue, where matrix elements are a predefined
#' quantile of the null model of given sig-sig interaction. If provided the
#' matrix will add another layer of filtering, where the absolute metric values
#' smaller than the provided value are set to 0.
#' @return A list with interaction metric matrices. The list element names match
#' those in the input network.lists.
#' @export

get_tissue_dataset_networks = function(tissue,
                                       network.lists,
                                       filter.list = NULL,
                                       filter.mat = NULL) {

    dump_zero_rows_cols = function(ll) {
      for (elem.name in names(ll) ) {
        mat = ll[[elem.name]]
        mat = mat[rowSums(abs(mat) ) > 0, colSums( abs (mat )) > 0]

        ll[[elem.name]] = mat
      }

      return(ll)
    }

    out.list = lapply(network.lists, function(x)
        x[[tissue]])

    if (is.null(filter.list)) {
        return( dump_zero_rows_cols (out.list) )
    }

    filter.list.absent = setdiff(names(filter.list), names(out.list))

    if (length( filter.list.absent > 0))  {
        stop(paste("filter.list has variable names not present in network.lists:",
                   paste(filter.list.absent, collapse = " ")) )
    }

    for (var in names(filter.list)) {
        varlim = filter.list[[var]]
        mat = out.list[[var]]
        mat[ abs(mat) < varlim ] = 0
        out.list[[var]] = mat
    }
    # For each matrix for metric thresholds, all the elements in the
    # calculated matrix that in absolute values are greater than the threshold
    # are set to 0
    if (! is.null(filter.mat)) {
      for (metric.name in names(filter.mat)) {
        out.list[[metric.name]] [ abs(out.list[[metric.name]]) < filter.mat[[metric.name]] ] = 0
      }
    }

    return( dump_zero_rows_cols (out.list) )
}
