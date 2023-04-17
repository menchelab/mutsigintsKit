#' Creates a tbl_graph network from a list of calculated metrics.
#'
#' @param network_list A list of calculated metrics matrices
#' @param network_names Names of metrics can be overriden with this parameter.
#'
#' @examples
#' mat = matrix(sample(100, 300, replace = TRUE), 30, 10,
#' 		   dimnames = list(paste0("c", 1:30), letters[1:10]))
#' mat[, 2] = 2 * mat[,1] + mat[,2]
#' mat[runif(300) < 0.5] = 0
#' 
#' 
#' metric_functions = list(CoDa = cor_coda,
#'                         cooccurrence = cooccurrence,
#'                         MI = bcmi,
#'                         Spearman = function(x) {cor_sigs(x, method = "spearman")})
#' 
#' metrics_out = lapply(metric_functions, function(x) x(mat))
#' graph_out = tissue_multi_graph(metrics_out)
#'
#' @import tidygraph
#'
#' @export
#' 
tissue_multi_graph = function(network_list, network_names = NULL) {
    
    if(! is.null(network_names)) {
        network_list = setNames(network_list, network_names)
    }
    
    for (netname in names(network_list)) {
        diag(network_list[[netname]]) = 0
        # excludes adjacency matrix where all elements are 0
        if ( ! sum( abs( network_list[[netname]] ) ) ) {
            network_list[[netname]] = NULL
        }
    }
    
    cat(names(network_list), "\n", sep = "\t")
    
    network_graph_list = lapply(names(network_list), function(net_name) {
        
        net = network_list[[net_name]]
        
        
        ### subject to being removed when as_tbl_graph bug is fixed
        if ( all(dim(net) == c(2,2) ) ) {
            net = cbind(rbind(net, 0), 0)
        }
        ### 
        
        tbl_net = as_tbl_graph(net, directed = FALSE) %>% 
            activate(edges) %>% 
            dplyr::mutate(method = net_name,
                   weight_ = weight,
                   weight = NULL)
        return(tbl_net)
    }) 
    
    network_graphs = Reduce(graph_join, network_graph_list)
    
    network_graphs = network_graphs %>% 
        activate(nodes) %>%
        filter(!node_is_isolated()) %>%
        activate(edges) %>% 
        dplyr::mutate(int_type = ifelse(weight_ > 0, "pos", "neg") ) %>% 
        dplyr::mutate(int_type = ifelse(method == "MI", "MI", int_type))
    return(network_graphs)
}
