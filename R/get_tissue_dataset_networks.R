#' Extracts the interaction metrics of a tissue 
#' from a list of interaction networks.
#' @param tissue Tissue to be extracted
#' @param network.lists The input list
#' @param filter.list A list of metric values to be filtered out. The list 
#' element names have the same names as metric names in network.lists. 
#' The values specify filtering threshold. Everything below the number is 
#' set to 0. Default: NULL
#' @return A list with interaction metric matrices. The list element names match
#' those in the input network.lists.
#' @export

get_tissue_dataset_networks = function(tissue,
                                       network.lists, 
                                       filter.list = NULL) {
    
    out.list = lapply(network.lists, function(x) 
        x[[tissue]])
    
    if (is.null(filter.list)) {
        return(out.list)
    } 
    
    filter.list.absent = setdiff(names(filter.list), names(out.list))
    
    if (length( filter.list.absent > 0))  {
        stop(paste("filter.llist has variable names not present:", 
                   paste(filter.list.absent, collapse = " ")) )
    } 
    
    for (var in names(filter.list)) {
        varlim = filter.list[[var]]
        mat = out.list[[var]]
        mat[ abs(mat) < varlim ] = 0
        mat = mat[rowSums(mat) > 0, colSums(mat) > 0]
        out.list[[var]] = mat
    }
    return(out.list)
}
