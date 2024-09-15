#' This function calculates interaction values using a given metric_func in a
#' dataset through bootstrapping.
#' @param data.input Input mutational signatures
#' @param metric_func Function calculating the interaction from a matrix
#' @param min.tissue.samples Tissues with less than
#' this number of samples are discarded. Default: 20
#' @param min.sig.mutcounts Signatures with total counts less than this value
#' will be filtered out. Default: 100
#' @param sample.rate The proportion of samples which should be sampled during
#' the bootstrapping procedure. Default: 0.9
#' @param sample.N Number of times to sample. Default: 100
#' @param N Number of times the function should be called for each sampled set:
#' Default: 1
#' @param seed Random number generator seed. Default: 1.
#' @param ... Parameters to be passed to metric_func.
#' @export

calculate_bootstrap_metric = function(data.input, metric_func,
                                      min.tissue.samples = 20,
                                      min.sig.mutcounts = 100,
                                      sample.rate = 0.9, sample.N = 100,
                                      N = 1, seed = 1,  ...) {

  action_func = function(x) {
    metric_func(x, ...)
  }

  tissues = unique(data.input$Cancer.Types)

  int.metric.list = list()

  for (tissue in tissues) {

    #### Filtering for tissues with enough samples
    tissue.signatures = list()
    cat("Processing tissue ", tissue, "\n")

    tissue.sig = data.input %>% filter(Cancer.Types == tissue) %>%
      dplyr::select(4:ncol(data.input))
    # tissue.sig = tissue.sig %>% filter(MMR == 0)
    tissue.sig = tissue.sig[, colSums(tissue.sig) > min.sig.mutcounts,
                            drop = FALSE]

    if (nrow(tissue.sig) < min.tissue.samples) {
      cat("Tissue:", tissue , " has <", min.tissue.samples,
          "samples in ", dataset, ". Skipping.\n")
      tissue.signatures = NULL
    } else {
      tissue.signatures = tissue.sig
    }

    metric_out = run_calc_metrics(tissue.signatures, action_func,
                                  sample.rate, sample.N, N, seed)
    int.metric.list[[tissue]] = metric_out
  }
  return(int.metric.list)
}
