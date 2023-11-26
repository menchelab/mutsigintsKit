#' The function aggregates mutations based on annotations
#' @param annotations annotations data frame. It contains must contain 2 columns-
#' Signature and Annotation
#' @param signature.matrix.df Signature matrix, which should be annotated
#' @return A matrix with the same number of rows as the input
#' siganture.matrix.df and aggregated and renamed columns
#' according to annotations dataframe.
#' @export


set_signature_annotation = function(annotations, signature.matrix.df) {

  sig.annotation.groups = annotations$Annotation [
    match(colnames(signature.matrix.df),
          annotations$Signature) ]

  for(i in 1:length(sig.annotation.groups)) {
    if (is.na(sig.annotation.groups[i]) ) {
      sig.annotation.groups[i] = colnames(signature.matrix.df)[i]
    }
  }

  signature.matrix.df.ann = signature.matrix.df %>% t() %>%
    rowsum(group = sig.annotation.groups) %>% t()

  return(signature.matrix.df.ann)
}

#' A function-wrapper to read sample sheet.
#' @param file.dir File directory.
#' @details The function can be later used in case some sort of processing
#' becomes necessary.
#'
#' @export

get_sample_sheet = function(file.dir) {
  read.delim(file.dir)
}


#' A function-wrapper to read specimen histology file.
#' @param file.dir File directory.
#' @details The function can be later used in case some sort of processing
#' becomes necessary.
#' @export

get_specimen_hist = function(file.dir) {
  read_excel(file.dir)
}


#' Performs p-value adjustment on a matrix of p-values based on input parameters.
#'
#' @return A matrix with the same dimensions and dimnames as the input
#'
p_matrix_process = function(p.values, p.adjust = TRUE, method = "BH") {
  if (p.adjust) {
    p.out <- p.values %>%
      as.vector %>% p.adjust(method = method) %>%
      matrix(ncol = ncol(p.values), dimnames = dimnames(p.values))
  } else {
    p.out = p.values
  }

}


#' Remove the columns and rows of the dataframe with all 0's
#' @param x input dataframe
#' @export
rm_zeros = function(x) {
  pos.x = abs(x)
  rmout = x[rowSums(pos.x) > 0, colSums(pos.x) > 0, drop = FALSE]
  return(rmout)
}


#' Significance starts from a pvalue.
#' @param p.val Input p value
get_sig_stars = function(p.val) {
  if (is.na(p.val)) {
    return(" ")
  }
  if (p.val < 0.001) {
    sig.star = "***"
  } else if (p.val < 0.01) {
    sig.star = "**"
  } else if (p.val < 0.05) {
    sig.star = "*"
  } else {
    sig.star = " "
  }
  return(sig.star)
}



#' Extract tissue and the signatures active in that tissue from a matrix of all
#' the data.
#'
#' @param sigs.full.subset is the signature activity matrix, where first column
#' is the cancer type/tissue and signature activities start from column 4. This
#' is the input type for signature activities from Alexandrov et al.(Nature, 2020)
#' paper.
#' @param tissue Tissue of interest
#'
#' @return Returns a subset of samples from the selected tissue and only active
#' signatures
#' @export

subset_tissue = function(sigs.full.subset, tissue) {

  tissue.signatures = sigs.full.subset %>% filter( Cancer.Types == tissue)

  signature.matrix = tissue.signatures[, 4:ncol(sigs.full.subset)]
  signature.matrix = signature.matrix[, which(colSums(signature.matrix) > 0),
                                      drop = FALSE]
  out = cbind(tissue.signatures[, 1:3], signature.matrix)
  return(out)
}


#' For a given signature activity matrix and pathway alteration matrix
#' computes Fisher's exact test for all pathway-signature interactions.
#' @param sigs.df signatures dataframe or matrix
#' @param pathways.df pathways dataframe or matrix
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A matrix with signature-number of rows and pathway-number of columns
#' with values indicating odds ratio of corresponding signatures and pathways.
#'
#' @export

get_sig_path_assocs = function(sigs.df, pathways.df, p.val.threshold = 0.05,
                               p.adjust = TRUE, method = "BH") {

  sigs = colnames(sigs.df)
  pathways = colnames(pathways.df)

  if (! all(rownames(sigs) == rownames(pathways))) {
    stop("Colnames corresponding to sample names should match between two matrices.")
  }

  odds.mat = matrix(0, nrow = length(sigs), ncol = length(pathways),
                    dimnames = list(sigs, pathways))

  p.values = matrix(0, nrow = length(sigs), ncol = length(pathways),
                    dimnames = list(sigs, pathways))

  for (sig in sigs) {
    for (pathway in pathways) {

      dd = cbind(sigs.df[, sig, drop = FALSE],
                 pathways.df[, pathway, drop = FALSE])
      dd = na.omit(dd)

      dd[dd > 0] = 1

      dd[, pathway] = factor(dd[ ,pathway], levels = c(0, 1),
                             labels = c("WT", "Mutated"))
      dd[ , sig] = dd[, sig] > 0
      dd[, sig] = factor(dd[, sig], levels = c(FALSE,TRUE),
                         labels = c("Inactive", "Active"))

      contingency.mat = table(dd)

      if (any(contingency.mat == 0)) {
        odds.mat[sig, pathway] = 0
        p.values[sig, pathway] = 1
      } else {
        fisher.out = fisher.test(contingency.mat) ### adding 1, to avoind +/- Inf values
        odds.mat[sig, pathway] = log2(fisher.out$estimate)
        p.values[sig, pathway] = fisher.out$p.value
      }
    }
  }

  p.out = p_matrix_process(p.values, p.adjust = p.adjust, method = method)

  # pheatmap(p.out)

  odds.mat [ p.out > p.val.threshold ] = 0

  # odds.mat = odds.mat[rowSums(abs(odds.mat), na.rm = TRUE) > 0, , drop = FALSE]
  # odds.mat = odds.mat[, colSums(abs(odds.mat), na.rm = TRUE) > 0 , drop = FALSE]

  return(odds.mat)
}



#' For signature and pathway matrices performs linear regression for
#' each signature-pathway pair.
#' For a given signature activity matrix and pathway alteration matrix
#' computes Fisher's exact test for all pathway-signature interactions.
#' @param sigs.df signatures dataframe or matrix
#' @param pathways.df pathways dataframe or matrix
#' @param sig.log Defines whether signatures should be logged. Default: TRUE
#' @param robust Specifies if robust statistics should be used.
#' Robust functions are called from robustbase package. Default: TRUE
#' @param path.to.sig If true, linear regression is used to predict signatures
#' from pathways. Otherwise, logistic regression is used to predict pathways
#' from signatures. Default: TRUE
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A matrix with signature-number of rows and pathway-number of columns
#' with values indicating odds ratio of corresponding signatures and pathways.
#' @export

get_sig_path_lms = function(sigs.df, pathways.df,
                            sig.log = TRUE,
                            robust = TRUE,
                            path.to.sig = TRUE,
                            p.val.threshold = 0.05,
                            p.adjust = TRUE, method = "BH",
                            ...) {

  if (sig.log) {
    sigs.df = sigs.df %>%
      mutate(across(.cols = everything(), ~ log(.x + 1 )))
  }

  sigs = colnames(sigs.df)
  pathways = colnames(pathways.df)

  if( ! all.equal(rownames(sigs), rownames(colnames)) ) {
    stop("The rownames (corresponding to sample names) in both datasets should match.")
  }

  tissue.concat = merge(pathways.df, sigs.df,
                        by = "row.names")

  int.mat = matrix(0, nrow = length(sigs), ncol = length(pathways),
                   dimnames = list(sigs, pathways))
  p.values = matrix(1, nrow = length(sigs), ncol = length(pathways),
                    dimnames = list(sigs, pathways))

  for (sig in sigs) {
    for (pathway in pathways) {
      # cat(sig, pathway, "\n")

      ## Applying a threshold on minimal number of non-zero elements
      # cat("before printing zero element counts.\n")
      zero.sigs = sum(tissue.concat[, sig,] != 0, na.rm = TRUE)
      zero.paths = sum(tissue.concat[, pathway] != 0, na.rm = TRUE)

      if (zero.sigs < 3 | zero.paths < 3) {
        # cat("zero.sigs = ", zero.sigs, " zero paths = ", zero.paths, "\n")
        int.mat[sig, pathway] = 0
        p.values[sig, pathway] = 1
        next
      }

      paths.binary = as.numeric(tissue.concat[, pathway] > 0)

      ### Drop if all elements are non-zero
      if (length(unique(paths.binary)) == 1) {
        int.mat[sig, pathway] = 0
        p.values[sig, pathway] = 1
        next
      }

      if ( robust ) {
        if (path.to.sig) {
          ## with Robustbase
          rob.lin.mod = robustbase::lmrob(
            tissue.concat[, sig] ~ tissue.concat[, pathway], ...)
          int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Estimate"][2]
          p.values[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Pr(>|t|)"][2]
          ## with rlm
          # rob.lin.mod = MASS::rlm(tissue.concat[, sig] ~ tissue.concat[, pathway], ...)
          # int.mat[sig, pathway] = summary(rob.lin.mod)$coefficients[, "Value"][2]
          # p.values[sig, pathway] = tryCatch({
          #     sfsmisc::f.robftest(rob.lin.mod, var = -1)$p.value},
          #     error = function(e) {return(1)})

        } else {
          ### Trying robust logistic regression first.
          ### If it throws an error, then regular logistic regression.

          log.mod = tryCatch({
            rob.out = robustbase::glmrob(
              paths.binary ~ 1 + tissue.concat[, sig],
              family = binomial)
            rob.out
          },
          error = function(e) {
            print(e$message)
            cat("Robust regression cannot be computed. \n Running regular regression instead.\n")
            glm.out = glm(
              paths.binary ~ 1 + tissue.concat[, sig],
              family = binomial)
            glm.out
          })
          int.mat[sig, pathway] = summary(log.mod)$coefficients[, "Estimate"][2]
          p.values[sig, pathway] = summary(log.mod)$coefficients[, "Pr(>|z|)"][2]
        }
      } else {
        if (path.to.sig) {
          # lin.mod = lm(tissue.concat[, sig] ~ tissue.concat[, pathway])
          lin.mod = lm(tissue.concat[, sig] ~ paths.binary)
          int.mat[sig, pathway] = summary(lin.mod)$coefficients[, "Estimate"][2]
          p.values[sig, pathway] = summary(lin.mod)$coefficients[,"Pr(>|t|)"][2]
        } else {
          paths.binary = as.numeric(tissue.concat[, pathway] > 0)
          log.mod = glm(
            paths.binary ~ 1 + tissue.concat[, sig],
            family = binomial, ...)
          int.mat[sig, pathway] = summary(log.mod)$coefficients[, "Estimate"][2]
          p.values[sig, pathway] = summary(log.mod)$coefficients[, "Pr(>|z|)"][2]
        }
      }
    }
  }

  p.out = p_matrix_process(p.values, p.adjust = p.adjust, method = method)

  int.mat[p.out > p.val.threshold] = 0
  # return(list(int.mat = int.mat, p.value = p.out))
  return(int.mat)
}

#' Extract signature and pathway activities for a tissue
#' @param tissue Tissue to extract
#' @param sigs.input Signature activities
#' @param pathways.input Pathway activities
#' @return A list with sigs and pathways elements with corresponding matrices.
#' Can be forwarded do interaction calculting functions
#' @export

get_tissue_pathway_activities = function(tissue,
                                         sigs.input,
                                         pathways.input) {

  tissue.sig.subset = subset_tissue(sigs.input, tissue = tissue)

  donor.ids = tissue.sig.subset$Sample.Names

  tissue.sig.subset = tissue.sig.subset %>% dplyr::select(4:ncol(.))

  # tissue.path.subset = pathways.input[ match(donor.ids,
  #                                            pathways.input$donor_id), ]

  tissue.path.subset = pathways.input %>% filter(donor_id %in% donor.ids) %>%
    arrange(factor(donor_id, levels = donor.ids) )

  tissue.path.subset = tissue.path.subset %>%
    dplyr::select(4:ncol(.)) %>%
    dplyr::select_if(colSums(., na.rm = TRUE) != 0)

  return(list(sigs = tissue.sig.subset, paths = tissue.path.subset))
}

#' Assess signature-pathway interactions across tissues for a custom function
#' @param sigs.input Signature activities. The rows correspond to samples, the
#' first three columns are Cancer.Types, Sample.Names, Accuracy, all the other
#' columns correspond to signature names. The values are signature activities.
#' @param pathways.input Pathway status - formatted like mutated.pathways.tissues
#' Pathway activities start at column 4. The first three columns are sample_id,
#' donor_id, Cancer.Types. Other columns correspond to pathway names. The values
#' correspond to number of mutations in the pathway.
#' @param interaction_function The function defining the metric. E.g. get_sig_path_assocs
#' @param path.min.tissues Minimal number of samples in each tissue to be considered
#' @param p.val.threshold p-value threshold for BH correction. Default: 0.05
#' @param p.adjust Controls if p-values should be adjusted. Default: TRUE
#' @param method P-value adjustement methods. Default: BH
#' @return A list with the length of abundant tissues in the datasets,
#' where each element is the interaction matrix
#' @export

sig_pathway_int = function(sigs.input,
                           pathways.input,
                           interaction_function,
                           path.min.tissues = 30,
                           p.val.threshold = 0.1,
                           p.adjust = TRUE,
                           method = "BH",
                           ...) {

  abundant.tissues = which(sigs.input$Cancer.Types %>%
                             table() > path.min.tissues) %>% names()

  tissue.odds.mats = list()

  for (tissue in abundant.tissues){
    cat(tissue, "\n")

    tissue.elems = get_tissue_pathway_activities(tissue,
                                                 sigs.input,
                                                 pathways.input)

    tissue.odds.mat = interaction_function(
      tissue.elems$sigs,
      tissue.elems$paths,
      p.val.threshold = p.val.threshold,
      p.adjust = p.adjust,
      method = method, ...) %>%
      rm_zeros()

    tissue.odds.mats[[tissue]] = tissue.odds.mat
  }
  return(tissue.odds.mats)
}


#' Assess signature-pathway interactions null model across tissues for a custom function
#' @param sigs.input Signature activities. The rows correspond to samples, the
#' first three columns are Cancer.Types, Sample.Names, Accuracy, all the other
#' columns correspond to signature names. The values are signature activities.
#' @param pathways.input Pathway status - formatted like mutated.pathways.tissues
#' Pathway activities start at column 4. The first three columns are sample_id,
#' donor_id, Cancer.Types. Other columns correspond to pathway names. The values
#' correspond to number of mutations in the pathway.
#' @param interaction_functions The functions defining the metrics.
#' E.g. get_sig_path_assocs
#' @param path.min.tissues Minimal number of samples in each tissue to be considered
#' @param N number of times the shuffling should take place
#' @return A list with the length of abundant tissues in the datasets,
#' where each element is a 3D array, where first two dims are sig and path, the 3rd
#' dim is the N metric values from N random shufflings.
#' @export

sig_pathway_int_null = function(sigs.input,
                                pathways.input,
                                interaction.functions,
                                path.min.tissues = 30,
                                N = 1000) {

  abundant.tissues = which(sigs.input$Cancer.Types %>%
                             table() > path.min.tissues) %>% names()

  tissue.odds.mats.null = list()

  shuffle_df = function(df) {
    out = data.frame(sapply(df, sample))
    rownames(out) = rownames(df)
    return(out)
  }

  for (tissue in abundant.tissues){
    cat(tissue, "\n")

    tissue.elems = get_tissue_pathway_activities(tissue,
                                                 sigs.input,
                                                 pathways.input)

    shuffled.sigs = replicate(N, shuffle_df(tissue.elems$sigs),
                              simplify = FALSE)
    ### no need to shuffled the pathways, since the sigs are already



    shuffled.metrics = lapply(
      interaction.functions,
      function(xf) {
        lapply(1:N, function(ind) {
          metric.out = xf(sigs = shuffled.sigs[[ind]],
                          paths = tissue.elems$paths)
          return(metric.out)
        })
      })

    metric.dists = list()
    for(metric.name in names(interaction.functions)) {
      metric.dists[[metric.name]] = simplify2array(shuffled.metrics[[metric.name]])
    }

    tissue.odds.mats.null[[tissue]] = metric.dists
  }
  return(tissue.odds.mats.null)
}


#' A ggheatmap wrapper for interaction matrix returns
#' @param metric.matrix the input matrix
#' @description All other parameters are passed to ggheatmap
#' \link[ggheatmap]{ggheatmap}
#' @return Returns a ggheatmap (which is a gg object,
#' but not a conventional gg object unfortunately.)
#' @export

ggheatmap_wrapper = function(metric.matrix,
                             colorPalette = "RdBu",
                             points = T, revColors = F, revRow = T,
                             scaleName = "log(OR)",
                             title = NULL,
                             limits = NULL,
                             ...) {

  if (( ! all(dim(metric.matrix) > 0) ) | is.null(dim(metric.matrix) ) ) {
    tissue.odds.plot = NA
    cat("The input matrix has no active dimensions or is null.")
  } else {
    orderrow = T
    ordercol = T
    if(dim(metric.matrix)[1] == 1) {
      orderrow = F
    }

    if(dim(metric.matrix)[2] == 1){
      ordercol = F
    }

    if (is.null(limits)) {
      range_lim = round(max(abs(range(metric.matrix))) * 1.1, 1)
      limits = c(-range_lim, range_lim)
    }

    tissue.odds.plot = myggheatmap(
      as.data.frame(metric.matrix), colorPalette = colorPalette,
      orderCol = ordercol,
      orderRow = orderrow,
      points = points,
      revColors = revColors,
      revRow = revRow,
      scaleName = scaleName,
      title = title,
      limits = limits, ...)
  }

  return(tissue.odds.plot)
}


#' #### TO BE REPLACED BY summarize_ints_to_df_clean
#' For a list of interaction metrics this function summarizes the positive and
#' negative interactions into a dataframe which can then be fed into plotting functions.
#'
#' @param list.of.int.elems A list with matrix elements.
#' @param threshold All the values below the threshold are discarded.
#' @param min.abssum Rows and columns containing with less than this number of
#' interactions (absolute) are removed. Default: 1.
#' @param upper.triangle.only If true, for symmetric matrices only upper triangle will be printed.
#' Default:TRUE.
#' @return Dataframe for plotting.
#' @export

summarize_ints_to_df = function(list.of.int.elems, threshold = 0.1, min.abssum = 1,
                                upper.triangle.only = TRUE) {

  all.rows = do.call(c, lapply(list.of.int.elems, function(x) rownames(x))) %>%
    unique()

  all.cols = do.call(c, lapply(list.of.int.elems, function(x) colnames(x))) %>%
    unique()

  pos.ints.mat = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                        dimnames = list(all.rows, all.cols))

  neg.ints.mat = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                        dimnames = list(all.rows, all.cols))

  for (rowelem in all.rows) {
    for (colelem in all.cols) {
      point.ints = unlist(sapply(list.of.int.elems,
                                 function(x) as.data.frame(x)[rowelem, colelem] ) )

      point.ints[abs(point.ints) < threshold] = 0

      pos.ints = sum(point.ints > 0, na.rm = TRUE)
      neg.ints = -sum(point.ints < 0, na.rm = TRUE)

      pos.ints.mat[rowelem, colelem] = pos.ints
      neg.ints.mat[rowelem, colelem] = neg.ints
    }
  }

  pos.pivot = pos.ints.mat %>%
    as.data.frame() %>% rownames_to_column(var = "rows") %>%
    pivot_longer(cols = -rows, values_to = "count", names_to = "cols") %>%
    mutate(int.type = "pos")

  # rows   name       count int.type
  # <chr>  <chr>      <dbl> <chr>
  # 1 Ageing Cell Cycle     1 pos
  # 2 Ageing HIPPO          2 pos
  # 3 Ageing NRF2           1 pos
  # 4 Ageing PI3K           2 pos

  neg.pivot = neg.ints.mat %>%
    as.data.frame() %>% rownames_to_column(var = "rows") %>%
    pivot_longer(cols = -rows, values_to = "count", names_to = "cols") %>%
    mutate(int.type = "neg")


  gg.final.dt = bind_rows(pos.pivot, neg.pivot)

  gg.final.dt = as.data.frame(gg.final.dt)

  abs.row.nonzero = gg.final.dt %>%
    group_by(rows) %>%
    summarise(abssum = sum(abs(count))) %>%
    filter(abssum >= min.abssum) %>%
    pull(rows)

  abs.col.nonzero = gg.final.dt %>%
    group_by(cols) %>%
    summarise(abssum = sum(abs(count))) %>%
    filter(abssum >= min.abssum) %>%
    pull(cols)


  row.indices = setNames(1:length(abs.row.nonzero), abs.row.nonzero)
  col.indices = setNames(1:length(abs.col.nonzero), abs.col.nonzero)

  gg.final.dt = gg.final.dt %>%
    filter(rows %in% abs.row.nonzero) %>%
    filter(cols %in% abs.col.nonzero) %>%
    mutate(x = row.indices[rows],
           y = col.indices[cols])

  ### If the matrices are symmetric, then only upper triangle is plotted.


  if (isTRUE(all.equal(neg.ints.mat, t(neg.ints.mat))) &
      isTRUE(all.equal(pos.ints.mat, t(pos.ints.mat))) & upper.triangle.only) {
    gg.final.dt = gg.final.dt %>% filter(x < y)

    row.indices = row.indices[ which(row.indices %in% gg.final.dt$x)]
    col.indices = col.indices[ which(col.indices %in% gg.final.dt$y)]
  }

  return(gg.final.dt)
}


#' For a list of interaction metrics this function summarizes the positive and
#' negative interactions, aggregates tissue names.
#'
#' @param int.value.list A list with matrix elements.
#' @param threshold All the values below the threshold are discarded.
#' @param min.abssum Rows and columns containing with less than this number of
#' interactions (absolute) are removed. Default: 1.
#' @return Dataframe for plotting.
#' @export

summarize_ints_to_df_new = function(int.value.list, threshold = 0.1, min.abssum = 1) {

  int.value.list = int.value.list[lengths(int.value.list) > 0]

  mat_to_df_list = lapply(names(int.value.list ),
                          function(nn) {
                            int.value.list[[nn]] %>%
                              as.data.frame() %>% rownames_to_column(var = "rows") %>%
                              pivot_longer(cols = -rows,
                                           values_to = "value", names_to = "cols") %>%
                              mutate(tissue = nn,
                                     int.type = ifelse(value > 0, "pos", "neg")) %>%
                              filter(value != 0, abs(value) > threshold)
                          })

  mat_to_df = do.call(rbind, mat_to_df_list)

  summarized = mat_to_df %>%
    group_by(rows, cols, int.type) %>%
    dplyr::summarize(count = n(), summedList = toString(sort(unique(tissue)))) %>%
    mutate(count = ifelse(int.type == "pos", count, count * (-1)))

  summarized = summarized %>% filter(abs(count) >= min.abssum)
  return(summarized)
}


#' For a list of interaction metrics this function summarizes the positive and
#' negative interactions in a plot.
#'
#' @param list.of.int.elems A list with matrix elements.
#' @param threshold All the values below the threshold are discarded.
#' @param min.abssum Rows and columns containing with less than this number of
#' interactions (absolute) are removed. Default: 1.
#' @param psize Controls the size of the triangles. Default: 8.
#' @param lsize Label size. Default: 2.
#' @param expand.mult A vector of expansion factors for x and y axes.
#' Default. c(0.04, 0.04)
#' @param upper.triangle.only If true, for symmetric matrices only diagonal will be printed.
#' Default:TRUE.
#'
#' @details The input is a list of matrices corresponding to e.g. interactions
#' in different tissues. The output plot summarizes the number of list elements
#' in which this interaction happens either in positive or negative direction.
#' The positive and negative interactions are summarized separately.
#'
#' @return ggplot object
#' @export


plot_all_counts = function(list.of.int.elems, threshold = 0.1, min.abssum = 1,
                           psize = 8, lsize = 2, expand.mult = c(0.04, 0.04),
                           upper.triangle.only = TRUE) {

  gg.final.dt = summarize_ints_to_df(list.of.int.elems, threshold, min.abssum)

  # gg.final.dt$text.col = sapply(gg.final.dt$count, function(x) {
  #     if ( x > 0 ) { return("black")
  #     }  else {return ("white")}
  # })
  # print(head(gg.final.dt))


  gg.final.dt = gg.final.dt %>%
    mutate(xlab = ifelse(int.type == "pos", x - 0.2, x + 0.2),
           ylab = ifelse(int.type == "pos", y + 0.15, y - 0.15),
           int.type = factor(int.type, levels = c("pos", "neg") ) )

  # gg.final.dt = gg.final.dt %>%
  #     mutate(xlab = x,
  #            ylab = y)

  gg.final.dt = as.data.frame(gg.final.dt)

  d <- ggplot(gg.final.dt, aes(x = x, y = y,
                               color = int.type,
                               shape = int.type,
                               label = abs(count),
                               alpha = abs(count)))

  if (Sys.info()['sysname'] == "Darwin") {
    # d = d + point_with_family(geom_point(size = 5), "wqy-microhei")
    d = d + point_with_family(geom_point(size = psize), "Arial Unicode MS")
  } else {
    d = d + geom_point(size = psize)
  }
  d = d +
    # point_with_family(geom_point(size = 18), "impact") +
    scale_shape_manual(values=c("\u25E4","\u25E2")) +
    scale_color_brewer(palette = "Set1") +
    geom_text(aes(x = xlab, y = ylab),
              size = lsize, fontface = "bold", color = "black")

  row.indices = gg.final.dt %>% select(rows, x) %>% unique()
  row.indices = setNames(row.indices$x, row.indices$rows) %>% base::sort()

  col.indices = gg.final.dt %>% select(cols, y) %>% unique()
  col.indices = setNames(col.indices$y, col.indices$cols) %>% base::sort()

  d = d +
    # theme_minimal() +
    theme( panel.border = element_blank(),
           panel.background = element_blank(),
           panel.grid = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           axis.ticks = element_blank(),
           # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
           axis.text.x = element_text(
             angle = 45,
             margin = margin(r = 0, t = 0, b = 0, l = 0),
             hjust = 0),
           axis.text.y = element_text(
             margin = margin(r = 0, t = 0, b = 0, l = 0)),
           axis.title = element_blank(),
           legend.position = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = expand.mult),
                       breaks = row.indices,
                       labels = names(row.indices),
                       position = "top") +
    scale_y_continuous(expand = expansion(mult = expand.mult),
                       breaks = col.indices,
                       labels = names(col.indices))
  d
  return(d)
}


# ### bipartite (needs more work)
# #' @param gg.df Dataframe returned by summarize_ints_to_df
# #' @return gggraph with sugiyama layout

# plot_bipartite = function(gg.df) {
#
#   df = summarize_ints_to_df(tissue.odds.ratio.unadjusted) %>% filter(count != 0)
#
#   # nodes.df
#
#   nodes.df = rbind(data.frame(id = unique(df$rows), type ="Sig"),
#                    data.frame(id = unique(df$cols), type = "Path"))
#
#   edges.df = df %>% rename(from = rows, to = cols)
#
#   bipart = tbl_graph(nodes = nodes.df, edges = edges.df)
#
#   pp = bipart %>% ggraph(layout = "sugiyama") +
#     geom_edge_parallel2(aes(edge_color = int.type, width = abs(count)),
#                         # width = 0.5,
#                         alpha = 0.8,
#                         start_cap = circle(0),
#                         end_cap = circle(0),
#                         sep = unit(0.8, 'mm')) +
#     geom_node_label(aes(label = id), fill = "white") +
#     theme_graph(base_family = "Arial") +
#     scale_edge_color_manual(values = c("pos" = "red4",
#                                        "neg" = "dodgerblue3") ) +
#     scale_edge_width_continuous(range = c(1,1.7)) #  +
#   # guides(width = "none")
#
#
#   dd = pp$data
#   dd$x = pp$data$y
#   dd$y = pp$data$x
#
#   pp$data = dd
#
#   # print(minor_plot(pp, expand.factor = 0.1))
#   return(pp)
# }


### bipartite (needs more work)
#' Bipartite network test
#' @param tissue.odds.ratios Dataframe returned by summarize_ints_to_df
#' @return gggraph with sugiyama layout
#' @export

plot_bipartite2 = function(tissue.odds.ratios) {

  df = summarize_ints_to_df_new(tissue.odds.ratios) %>% filter(count != 0)

  # nodes.df

  nodes.df = rbind(data.frame(id = unique(df$rows), type ="Sig"),
                   data.frame(id = unique(df$cols), type = "Path"))

  edges.df = df %>% rename(from = rows, to = cols)

  bipart = tbl_graph(nodes = nodes.df, edges = edges.df)

  bipart.nodes = bipart %>%
    activate(nodes) %>% pull(id) %>% unique()

  bipart = bipart %>%
    activate(nodes) %>%
    mutate(annot.class = ifelse(type == "Path", "Pathway", signature.annotations %>%
                                  filter(Annotation %in% bipart.nodes) %>%
                                  select(Origin, Annotation) %>%
                                  unique() %>%
                                  arrange(factor(Annotation, levels = bipart.nodes)) %>%
                                  pull(Origin)) )

  calc.layout = create_layout(bipart, layout = "sugiyama")

  pp = ggraph(graph = calc.layout) +
    geom_edge_parallel(aes(edge_color = int.type, width = abs(count)),
                       # width = 0.5,
                       alpha = 0.7,
                       start_cap = circle(0),
                       end_cap = circle(0),
                       sep = unit(0.8, 'mm')) +
    geom_node_point(aes(fill = annot.class, shape = type, color = annot.class),
                    # color = "gray30",
                    size = 4,
                    stroke = 0.6) +
    scale_shape_manual(values = c("Path" = 22,
                                  "Sig" = 21),
                       labels = c("Pathway", "Signature")) +
    scale_color_manual(values = c("Endogenous" = "#CAE7B9",
                                  "Environmental" = "#f3de8a",
                                  "Clock-like" = "white",
                                  "Unknown" = "#4e6151",
                                  "Pathway" = "gray50"),
                       name = "Signature class") +
    scale_fill_manual(values = c("Endogenous" = "#CAE7B9",
                                 "Environmental" = "#f3de8a",
                                 "Clock-like" = "white",
                                 "Unknown" = "#4e6151",
                                 "Pathway" = "gray50"),
                      name = "Signature class") +
    geom_node_point(aes(shape = type),
                    color = "gray30",
                    size = 4,
                    stroke = 0.6) +
    geom_node_label(aes(label = id#, fill = type
    ),
    color = "gray10",
    label.size = 0,
    hjust = ifelse(calc.layout[,2] > 1, -0.45, 1.3),
    label.r = unit(0.2, "lines"),
    label.padding = unit(0, "lines"),
    nudge_y = 0) +
    theme_graph(base_family = "Arial") +
    scale_edge_color_manual(values = c("pos" = "red4",
                                       "neg" = "dodgerblue3"),
                            labels = c("positive", "negative")) +
    scale_edge_width_continuous(breaks = seq(bipart %>%
                                               activate(edges) %>%
                                               pull(count) %>%
                                               abs %>%max),
                                range = c(0.8,1.5),
                                # guide = "none"
                                name = "Count") #  +
  # guides(width = "none")

  dd = pp$data
  dd$x = pp$data$y
  dd$y = pp$data$x

  pp$data = dd

  # print(minor_plot(pp, expand.factor = 0.1))
  return(pp)
}




#' Summarizing individual matrices of interactions
#'
#' @param interaction.list Input list
#' @param pos.ints Defines if positive or negative interactions
#' should be considered. Default: TRUE
#' @return Matrix with element is concatenated tissue names
#' bearing the interaction
#'
summarize_int_mat = function(interaction.list, pos.ints = TRUE) {

  interaction.list = lapply(interaction.list, as.data.frame)

  active.colnames = do.call(c, sapply(interaction.list, colnames) ) %>%
    unique() %>% sort
  active.rownames = do.call(c, sapply(interaction.list, rownames) ) %>%
    unique() %>% sort

  outmat = matrix("", nrow = length(active.rownames),
                  ncol = length(active.colnames),
                  dimnames = list(active.rownames,
                                  active.colnames))

  for(rowelem in active.rownames) {
    for(colelem in active.colnames) {
      vals = sapply(interaction.list, function(x)
        as.data.frame(x)[rowelem, colelem]) %>% unlist
      if(pos.ints) {
        mat.names = which(vals > 0) %>% names
      } else {
        mat.names = which(vals < 0) %>% names
      }
      outmat [rowelem, colelem] = paste(mat.names, collapse = ",")
    }
  }
  return(outmat)
}

### https://stackoverflow.com/questions/48531257/ggplot-geom-point-how-to-set-font-of-custom-plotting-symbols
point_with_family <- function(layer, family) {
  old_geom <- layer$geom
  new_geom <- ggproto(
    NULL, old_geom,
    draw_panel = function(self, data, panel_params, coord, na.rm = FALSE) {
      pts <- ggproto_parent(GeomPoint, self)$draw_panel(
        data, panel_params, coord, na.rm = na.rm
      )
      pts$gp$fontfamily <- family
      pts
    },
    draw_key = function(self, data, params, size) {
      pts <- ggproto_parent(GeomPoint, self)$draw_key(
        data, params, size
      )
      pts$gp$fontfamily <- family
      pts
    }
  )
  layer$geom <- new_geom
  layer
}


#' Calculates interaction values using a given metric_func in a dataset
#' through bootstrapping.
#' @param data.input Input mutational signatures
#' @param metric_func Function calculating the interaction from a matrix
#' @param min.tissue.samples Tissues with less than
#' this number of samples are discarded. Default: 20
#' @param sample.rate The proportion of samples which should be sampled during
#' the bootstrapping procedure. Default: 0.9
#' @param sample.N Number of times to sample. Default: 100
#' @param N Number of times the function should be called for each sampled set:
#' Default: 1
#' @param seed Random number generator seed. Default: 1.
#' @param ... Parameters to be passed to metric_func.
#' @export

get_metrics_list = function(data.input, metric_func,
                            min.tissue.samples = 20,
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

    dt = data.input
    tissue.sig = dt %>% filter(Cancer.Types == tissue) %>%
      dplyr::select(4:ncol(dt))
    # tissue.sig = tissue.sig %>% filter(MMR == 0)
    tissue.sig = tissue.sig[, colSums(tissue.sig) > 10,
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


# Calculates the metric for a list or a dataframe.

run_calc_metrics = function(tissue.signatures, action_func, sample.rate = 0.9,
                            sample.N = 100, N = 1, seed = 1) {

  if (inherits(tissue.signatures, "list")) {
    int.metric.calcs = lapply(tissue.signatures, function(dt)
      calc_metric(dt, action_func, sample.rate, sample.N, N, seed)
    )
  } else if (inherits(tissue.signatures, "data.frame")) {
    int.metric.calcs = calc_metric(tissue.signatures, action_func,
                                   sample.rate, sample.N, N, seed)
  }
  return(int.metric.calcs)
}

#' Calculates the metric function for a matrix, by performing bootstrapping.
#' @param dt.sigs A dataframe or a matrix of mutational signatures.
#' @param action_func Function calculating the interaction from a matrix
#' @param sample.rate The proportion of samples which should be sampled during
#' the bootstrapping procedure. Default: 0.9
#' @param sample.N Number of times to sample. Default: 100
#' @param N Number of times the function should be called for each sampled set:
#' Default: 1
#' @param seed Random number generator seed. Default: NULL.

calc_metric = function(dt.sigs, action_func, sample.rate, sample.N, N, seed = NULL ) {

  if (!is.null(seed)) { set.seed(seed)}
  if (is.null(dt.sigs)) {
    return(NULL)
  } else {
    sample.length = nrow(dt.sigs)
    sample.counts = round(sample.length * sample.rate)

    sampled.out = lapply(1:sample.N, function(i) {
      sampled.signatures = sample_n(dt.sigs, sample.counts)
      ll <- replicate(N, action_func(sampled.signatures),
                      simplify = FALSE )

      ll.named = lapply(ll, function(x) {
        colnames(x) = colnames(sampled.signatures)
        rownames(x) = colnames(sampled.signatures)
        return(x)
      })
      return(ll.named)

    }) }
}



#' Plotting a heatmap for signature networks. Can be a PCAWG-format with first
#' 3 columns being a metadata, or a matrix/data.frame with all signature activity
#' values.
#' @param .dataset Input dataset
#' @param filename If provided the pheatmap will be save with this filename.
#' Default: NULL
#' @param main Title of the heatmap
#' @param border_color Border color of cells. This parameter is there to control
#' removing border color with NA. Pheatmap doesn't properly remove it
#' and the border is still present when saving the plot. Default: gray60.
#' Should be NA to remove the border_color.
#' @param ... Parameters will be passed to pheatmap function.
#'
#' @export

plot_tissue_heatmap = function(.dataset, filename = NULL,
                               main = NULL, border_color = "grey60", ...) {

  classes = sapply(.dataset, class)
  if (! all(sapply(.dataset, class)[1:3] %in% c("numeric", "integer")))  {
    dt.plot = .dataset[, 4:ncol(.dataset)]
  } else {

    dt.plot = .dataset
  }

  dt.plot = dt.plot[, colSums(dt.plot) > 0]

  if (is.null(main)) {
    main = paste(nrow(.dataset), "samples")
  } else {
    main = paste(main, "-", nrow(.dataset), "samples")
  }

  if (is.null(filename) ) {
    pp = pheatmap(log(dt.plot + 1),
                  color = viridis(15),
                  show_rownames = FALSE,
                  main = main,
                  width = 7.1, height = 6.83,
                  border_color = border_color,
                  ...)
  } else {
    pp = pheatmap(log(dt.plot + 1),
                  color = viridis(15),
                  show_rownames = FALSE,
                  width = 7.1, height = 6.83,
                  main = main,
                  filename = filename,
                  border_color = border_color,
                  ...)
  }

  if (is.na(border_color)) {
    # https://stackoverflow.com/questions/44318690/no-border-color-in-pheatmap-in-r
    grob_classes <- purrr::map(pp$gtable$grobs, class)
    idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(pp$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]

    ## Remove borders around cells
    pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
    pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0
  }
  return(pp)
}


#' Adds a legend title to pheatmap legend
#'
#' @param p ggplot object produced by pheatmap.
#' @param legend.grob.index The grob which contains the legend grob.
#' @details Based on Mike H.'s answer from https://stackoverflow.com/questions/36852101/r-legend-title-or-units-when-using-pheatmap
#' @export

add_pheatmap_legend_title = function(p, legend.text = "log(n)") {

  legend.grob.index = which(p$gtable$layout$name == "legend")

  legend.grob <- p$gtable$grob[[legend.grob.index]]

  legend.grob$children[[1]]$y <- legend.grob$children[[1]]$y - 0.08 * legend.grob$children[[1]]$y
  legend.grob$children[[2]]$y <- legend.grob$children[[2]]$y - 0.08 * legend.grob$children[[2]]$y
  legend.grob$children[[1]]$x <- legend.grob$children[[1]]$x + 0.15 * legend.grob$children[[1]]$x
  legend.grob$children[[2]]$x <- legend.grob$children[[2]]$x + 0.15 * legend.grob$children[[2]]$x

  leg_label <- textGrob(legend.text,
                        x=0,y=0.95,hjust=0,vjust=0,gp=gpar(fontsize=10,fontface="bold"))

  legend.grob2 <- addGrob(legend.grob,leg_label)

  p$gtable$grobs[[legend.grob.index]] = legend.grob2
  return(p)
}


#' The function summarizes interactions for a given metric across tissues,
#' identifies common interaction motifs.
#' @param metric.list input metric.list, element names are different interaction
#' metrics. Each element is a list corresponding to interactions corresponding
#' to individual tissues.
#' @param metric Metric for which the summaries should be obtained.
#' @param outdir The output directory where the interaction matrices and the
#' outputs should be written.
#' @param threshold All interactions below this threshold are set to 0.
#' Default: 0.2.
#' @export

get_common_sigs = function(metric.list, metric, outdir, threshold = 0.2) {

  for (tissue in names(metric.list[[metric]])) {
    cat("\t\tTissue = ", tissue, "\n")
    tissue.int.mat = metric.list[[metric]][[tissue]]

    if(is.null(tissue.int.mat)) next

    tissue.int.mat[ abs(tissue.int.mat) < threshold] = 0
    tissue.int.mat  = sign(tissue.int.mat)
    write.table(tissue.int.mat, file = paste0(outdir, "/",  tissue, ".txt"),
                sep = "\t", quote = F, row.names = TRUE, col.names = TRUE)
  }
  python.call = paste("python3", here("python", "combine_networks.py"),
                      outdir)
  cat(python.call, "\n")
  system(python.call, wait = TRUE)
}


#' Summarize different metric outputs for an interaction
#' @param network.lists The input list.
#' @param tissue Tissue to be extracted.
#' @param filter.list A list of metric values to be filtered out. The list
#' element names have the same names as metric names in network.lists.
#' The values specify filtering threshold. Everything below the number is
#' set to 0. Default: NULL
#' @param filter.mat Passed to get_tissue_dataset_networks.
#' @export

concat_networks = function(network.lists, tissue, filter.list, filter.mat = NULL) {

  try({
    tissue.nets = get_tissue_dataset_networks(
      tissue,
      network.lists = network.lists,
      filter.list = filter.list,
      filter.mat = filter.mat
    )

    for (type in names(tissue.nets)) {
      if (sum(abs(tissue.nets[[type]]) ) == 0  ) {
        tissue.nets[[type]] = NULL
      }
    }

    ts.multi.graph = tissue_multi_graph(tissue.nets)
    ## layout = "dh"
    ts.adj.mat = as.matrix(as_adjacency_matrix(ts.multi.graph))
    ts.adj.mat[ts.adj.mat < 2] = 0
    ts.adj.mat[ts.adj.mat > 0] = 1
    ts.adj.mat = ts.adj.mat + t(ts.adj.mat)


    edge.info = ts.multi.graph %>% activate(edges) %>% data.frame() %>%
      filter(int_type != "MI")

    adj.net = as.matrix(as_adjacency_matrix(ts.multi.graph)) * 0
    for(i in 1:nrow(edge.info)) {
      if (adj.net[edge.info[i, "from"], edge.info[i, "to"]] != 0 ) {
        if (adj.net[edge.info[i, "from"], edge.info[i, "to"]] *  edge.info[i, "weight_"] < 0) {
          cat("\n####################i = ",i, " the signs differ.\n")
        }
      }
      adj.net[edge.info[i, "from"], edge.info[i, "to"]] =
        adj.net[edge.info[i, "from"], edge.info[i, "to"]] + edge.info[i, "weight_"]
    }
    adj.net = adj.net + t(adj.net)
    ts.adj.mat = ts.adj.mat * adj.net
    ts.adj.mat = ts.adj.mat [colSums(abs(ts.adj.mat) ) > 0,
                             rowSums(abs(ts.adj.mat)) > 0]
    return(ts.adj.mat)
  })

  return(c())
}



#' Wrapper of tissue interaction network plot
#' @param tissue Tissue to be extracted
#' @param network.lists The input list
#' @param filter.list A list of metric values to be filtered out. The list
#' element names have the same names as metric names in network.lists.
#' The values specify filtering threshold. Everything below the number is
#' set to 0. Default: NULL
#' @param filter.mat Passed to the function get_tissue_dataset_networks
#' @param layout layout to be passed to ggraph. Default: stress. Some useful
#' alternative is dh.
#' @export

plot_multi_network = function(network.lists, tissue,
                              filter.list = NULL,
                              filter.mat = NULL, layout = "stress") {

  tissue.nets = get_tissue_dataset_networks(
    tissue,
    network.lists = network.lists,
    filter.list = filter.list,
    filter.mat = filter.mat
  )

  for (type in names(tissue.nets)) {
    if (sum(abs(tissue.nets[[type]]) ) == 0  ) {
      tissue.nets[[type]] = NULL
    }
  }

  ts.multi.graph = tissue_multi_graph(tissue.nets)
  pp = print_multi_graphs(ts.multi.graph, layout = layout)

}


#' Summarize positive and negative interactions across tissues
#'
#' For each signature pair the function counts how many positive and negative
#' interactions were observed across tissues and returns a matrix of lists with
#' two elements called pos and neg or NULL if no interactions were found.
#'
#' @param summary.list A list of interactions for individual tissues. Each
#' element is a matrix.
#' @param tissue.names Input tissue names
#' @return a list of two matrices for each dataset. The elements of the matrices
#' are lists as elements with two elements called pos and neg for positive and
#' negative interactions.

summarize_by_tissues = function(summary.list, tissue.names = FALSE) {

  all.signatures = sapply(networks.concat, rownames) %>%
    unlist() %>%
    unique()

  all.mat = matrix(list(), ncol = length(all.signatures),
                   nrow = length(all.signatures),
                   dimnames = list(all.signatures, all.signatures))

  for(tissue in names(summary.list)) {
    cat("tissue =", tissue, "\n")
    dt_res = summary.list[[tissue]]

    signames = colnames(dt_res)

    if (is.null(dt_res)) {
      dt.inds.pos = c()
      dt_inds_neg = c()
      next
    }

    dt.inds.pos = which(dt_res > 0, arr.ind = TRUE)

    if ( length(dt.inds.pos) != 0) {
      for(i in 1:nrow(dt.inds.pos)) {
        sig1 = signames[dt.inds.pos[i, "row"]]
        sig2 = signames[dt.inds.pos[i, "col"]]
        if (tissue.names) {
          if (is.null(all.mat[sig1, sig2][[1]]))  {
            all.mat[sig1, sig2][[1]] = list(pos = "", neg = "")
          }
          all.mat[sig1, sig2][[1]][["pos"]] =
            paste(tissue, all.mat[sig1, sig2][[1]][["pos"]], sep = ",")

        } else {
          if (is.null(all.mat[sig1, sig2][[1]]))  {
            all.mat[sig1, sig2][[1]] = list(pos = 0, neg = 0)
          }
          all.mat[sig1, sig2][[1]][["pos"]] = all.mat[sig1, sig2][[1]][["pos"]]  + 1
        }
      }
    }

    dt.inds.neg = which(dt_res < 0, arr.ind = TRUE)
    if(length(dt.inds.neg) != 0 ) {
      for(i in 1:nrow(dt.inds.neg)) {
        sig1 = signames[dt.inds.neg[i, "row"]]
        sig2 = signames[dt.inds.neg[i, "col"]]
        if (tissue.names) {
          if (is.null(all.mat[sig1, sig2][[1]])){
            all.mat[sig1, sig2][[1]] = list(pos = "", neg = "")
          }
          all.mat[sig1, sig2][[1]][["neg"]] =
            paste(tissue, all.mat[sig1, sig2][[1]][["neg"]], sep = ",")
        } else {
          if (is.null(all.mat[sig1, sig2][[1]])){
            all.mat[sig1, sig2][[1]] = list(pos = 0, neg = 0)
          }
          all.mat[sig1, sig2][[1]][["neg"]] = all.mat[sig1, sig2][[1]][["neg"]] - 1
        }
      }
    }
  }

  return(all.mat)
}


#' For a list of interaction metrics this function summarizes the positive and
#' negative interactions in a matrix, where for each interaction the cell
#' contains a comma-separated tissue  names where this interaction was observed.
#'
#' @param list.of.int.elems A list with matrix elements.
#' @param threshold All the values below the threshold are discarded.
#'
#' @details The input is a list of matrices corresponding to e.g. interactions
#' in different tissues. The output plot summarizes the number of list elements
#' in which this interaction happens either in positive or negative direction.
#' The positive and negative interactions are summarized separately.
#'
#' @return ggplot object
#' @export

get_interaction_tissues = function(list.of.int.elems, threshold = 0.1) {

  # all.sigs = do.call(c, lapply(summary.matrix, function(x) colnames(x)) ) %>%
  #     unique()

  all.rows = do.call(c, lapply(list.of.int.elems, function(x) rownames(x))) %>%
    unique()

  all.cols = do.call(c, lapply(list.of.int.elems, function(x) colnames(x))) %>%
    unique()

  pos.ints.tissues = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                            dimnames = list(all.rows, all.cols))

  neg.ints.tissues = matrix(0, nrow = length(all.rows), ncol = length(all.cols),
                            dimnames = list(all.rows, all.cols))

  for (rowelem in all.rows) {
    for (colelem in all.cols) {
      point.ints = unlist(sapply(list.of.int.elems,
                                 function(x) as.data.frame(x)[rowelem, colelem]))

      point.ints[abs(point.ints) < threshold] = 0

      pos.ints = sum(point.ints > 0, na.rm = TRUE)
      neg.ints = -sum(point.ints < 0, na.rm = TRUE)

      pos.tissues = paste(names(which(point.ints > 0)), collapse = ", ")
      neg.tissues = paste(names(which(point.ints < 0)), collapse = ", ")


      pos.ints.tissues[rowelem, colelem] = pos.tissues
      neg.ints.tissues[rowelem, colelem] = neg.tissues
    }
  }

  pos.ints.tissues = pos.ints.tissues[rowSums(pos.ints.tissues > "") > 0,
                                      colSums(pos.ints.tissues > "") > 0 ]

  neg.ints.tissues = neg.ints.tissues[rowSums(neg.ints.tissues > "") > 0,
                                      colSums(neg.ints.tissues > "") > 0 ]

  return(list(pos.tissues = pos.ints.tissues,
              neg.tissues = neg.ints.tissues))
}


#' Get survival data for relevant tissues
#' @param dataset Signature-activities for samples. E.g.PCAWG.full.subset.ann.
#' The first three columns are Cancer.Types, Sample.Names, Accuracy and fourth on
#' signature activities.
#' @param clin.df Clinical data for samples. E.g. PCAWG.clin.df.
#' @param tissues A vector of tissues where the interaction of signatures should
#' be assessed.
#' @param is.TCGA Indicate if TCGA data is used
#' @export
get_relevant_clin_df = function(clin.df, dataset, is.TCGA, tissues) {

  tissues.subset = dataset %>% filter(Cancer.Types %in% tissues)

  tissues.subset = tissues.subset[ ! duplicated(tissues.subset$Sample.Names), ]

  if ( ! is.TCGA ) { ### for PCAWG
    relevant.clin.df = clin.df[ match(tissues.subset$Sample.Names,
                                      clin.df$icgc_donor_id), ]

    non_na.indices = which(!is.na(relevant.clin.df$icgc_specimen_id))

    tissues.subset = tissues.subset[non_na.indices, ]
    relevant.clin.df = relevant.clin.df[ non_na.indices, ]

  } else { ### for TCGA

    tissues.subset$Sample.Names = substr(tissues.subset$Sample.Names, 1, 12)
    valid.indices = which(tissues.subset$Sample.Names %in% clin.df$bcr_patient_barcode)
    tissues.subset = tissues.subset[valid.indices, ]

    relevant.clin.df = clin.df[ match(tissues.subset$Sample.Names,
                                      clin.df$bcr_patient_barcode), ]

    relevant.clin.df = as_tibble(relevant.clin.df)

    relevant.clin.df = rename(relevant.clin.df, survival_time = times)
    relevant.clin.df = rename(relevant.clin.df, vital_status = patient.vital_status)
  }

  return(list(relevant.clin.df = relevant.clin.df, tissues.subset = tissues.subset))
}


#' Survival analysis for signature-signature interactions with Cox proportional
#' hazards regression model.
#'
#' @param dataset Signature-activities for samples. E.g.PCAWG.full.subset.ann.
#' The first three columns are Cancer.Types, Sample.Names, Accuracy and fourth on
#' signature activities.
#' @param clin.df Clinical data for samples. E.g. PCAWG.clin.df.
#' @param signatures A vector of two signatures which should be considered for
#' survival analysis.
#' @param tissues A vector of tissues where the interaction of signatures should
#' be assessed.
#' @param legend_pos Legend position. Default: c(0.8, 0.8).
#' @param with.total.muts If TRUE the total number of mutations in the samples will
#' be provided as a confounder to the model. Default: TRUE
#' @param tmb.logged If TRUE the tumor mutational burden will be logged.
#' Default: TRUE
#' @param binary.status If TRUE, the model will compare samples with both signatures
#' with all the other samples having either of the signatures or none. Default:
#' FALSE
#' @param epistatic The survival model tests the interaction of int1 and int2 as
#' int1 * int2. Default: FALSE
#' @param get_df If true, the function returns the constructed model dataframe
#' (for debugging). Default:FALSE.
#' @param conf.int If TRUE confidence interval should be added to ggsurvplot.
#' @export

survival_for_interactions = function(dataset, clin.df, signatures,
                                     tissues, legend_pos = c(0.8, 0.8),
                                     with.total.muts = TRUE,
                                     tmb.logged = TRUE,
                                     binary.status = FALSE,
                                     epistatic = FALSE,
                                     get_df = FALSE,
                                     conf.int = FALSE) {

  # SBS40_APOBEC = survival_for_interactions(dataset = PCAWG.full.subset.ann,
  #                                          signatures = c("SBS40", "APOBEC"),
  #                                          tissues = c("Ovary-AdenoCA", "Thy-AdenoCA"),
  #                                          clin.df = PCAWG.clin.df)


  #### checking if TCGA
  TCGA = FALSE
  if ( substr(dataset$Sample.Names[1], 1, 4) == "TCGA") {
    TCGA = TRUE
  }

  if (length(signatures) != 2) {
    stop(paste(c("2 signatures are required for this function:", signatures) ),
         collapse = " ")
    return(NULL)
  }

  relevant.clin.df.out = get_relevant_clin_df(clin.df = clin.df, dataset = dataset,
                                              is.TCGA = TCGA, tissues = tissues)

  relevant.clin.df =  relevant.clin.df.out$relevant.clin.df
  tissues.subset = relevant.clin.df.out$tissues.subset
  # tissues.subset = dataset %>% filter(Cancer.Types %in% tissues)
  #
  # tissues.subset = tissues.subset[ ! duplicated(tissues.subset$Sample.Names), ]
  #
  # if ( ! TCGA ) { ### for PCAWG
  #     relevant.clin.df = clin.df[ match(tissues.subset$Sample.Names,
  #                                       clin.df$icgc_donor_id), ]
  #
  #     non_na.indices = which(!is.na(relevant.clin.df$icgc_specimen_id))
  #
  #     tissues.subset = tissues.subset[non_na.indices, ]
  #     relevant.clin.df = relevant.clin.df[ non_na.indices, ]
  #
  # } else { ### for TCGA
  #
  #     tissues.subset$Sample.Names = substr(tissues.subset$Sample.Names, 1, 12)
  #     valid.indices = which(tissues.subset$Sample.Names %in% clin.df$bcr_patient_barcode)
  #     tissues.subset = tissues.subset[valid.indices, ]
  #
  #     relevant.clin.df = clin.df[ match(tissues.subset$Sample.Names,
  #                                       clin.df$bcr_patient_barcode), ]
  #
  #     relevant.clin.df = as_tibble(relevant.clin.df)
  #
  #     relevant.clin.df = rename(relevant.clin.df, survival_time = times)
  #     relevant.clin.df = rename(relevant.clin.df, vital_status = patient.vital_status)
  #
  # }

  survival.df = cbind(tissues.subset[, 1:3], tissues.subset[, signatures],
                      relevant.clin.df[, c("survival_time", "vital_status",
                                           "age_at_diagnosis"#, "sex"
                      ) ] )

  survival.df = cbind(survival.df, total_muts =
                        rowSums(tissues.subset[, 4:ncol(tissues.subset)]))

  signatures = sort(signatures)
  sig1 = signatures[1]
  sig2 = signatures[2]

  colnames(survival.df)[ which(colnames(survival.df) == sig1)] = "SBS__1"
  colnames(survival.df)[ which(colnames(survival.df) == sig2)] = "SBS__2"

  survival.df = survival.df %>% filter(!is.na(SBS__1), !is.na(SBS__2))

  survival.df$exists__1 = as.numeric(survival.df$SBS__1 > 0)
  survival.df$exists__2 = as.numeric(survival.df$SBS__2 > 0)
  survival.df$exists__12 = as.numeric(survival.df$SBS__1 > 0 & survival.df$SBS__2 > 0)
  survival.df$exists__None = as.numeric(survival.df$SBS__1 == 0 & survival.df$SBS__2 == 0)


  sig.comb.status = apply(survival.df[, c("exists__1", "exists__2")], MARGIN = 1,
                          function(x) {
                            if(x[1] == 1 & x[2] == 1) return(paste0(sig1, "+", sig2))
                            if(x[1] == 0 & x[2] == 1) return(paste0( sig2))
                            if(x[1] == 1 & x[2] == 0) return(paste0( sig1))
                            if(x[1] == 0 & x[2] == 0) return(paste0("None"))
                          } )

  survival.df$status = sig.comb.status
  survival.df$status = factor(survival.df$status, levels = c("None", sig1,
                                                             sig2, paste0(sig1, "+", sig2)))


  if (binary.status) {
    survival.df$status = ifelse(survival.df$status == paste0(sig1, "+", sig2),
                                paste0(sig1, "+", sig2), "Others")
    survival.df$status = factor(survival.df$status,
                                levels = c("Others", paste0(sig1, "+", sig2)))

  }
  if(get_df) {
    return(survival.df)
  }


  formula.string = "Surv(survival_time, vital_status) ~ age_at_diagnosis"

  if (length(tissues) > 1) {
    formula.string = paste0(formula.string, " + Cancer.Types")
  }

  if (with.total.muts) {
    if (tmb.logged) {
      formula.string = paste0(formula.string, " + log(total_muts + 1)")
    } else {
      formula.string = paste0(formula.string, " + total_muts")
    }
  }

  if (epistatic) {
    formula.string = paste0(formula.string, " + exists__1 + exists__2 + exists__1 * exists__2")
  } else {
    formula.string = paste0(formula.string, " + status")
  }

  cox <- coxph(as.formula(formula.string), data = survival.df, na.action = na.omit,
               x = TRUE)

  if (epistatic) {
    nn = names(cox$coefficients)
    nn[nn == "exists__1"] = sig1
    nn[nn == "exists__2"] = sig2
    nn[nn == "exists__1:exists__2"] = paste(sig1, sig2, sep = "*")
    names(cox$coefficients) = nn
  }
  # if (length(tissues) >1 ) {
  #   if (with.total.muts) {
  #     if (tmb.logged) {
  #       cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                      Cancer.Types + status + log(total_muts + 1),
  #                    data = survival.df, na.action = na.omit)
  #     } else {
  #       cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                      Cancer.Types + status + total_muts,
  #                    data = survival.df, na.action = na.omit)
  #     }
  #   }else {
  #     cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                    Cancer.Types + status,
  #                  data = survival.df, na.action = na.omit)
  #   }
  # } else {
  #   if (with.total.muts) {
  #     if (tmb.logged) {
  #       cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                      status + log(total_muts + 1),
  #                    data = survival.df, na.action = na.omit)
  #     } else {
  #       cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                      status + total_muts,
  #                    data = survival.df, na.action = na.omit)
  #
  #     }
  #   }else {
  #     cox <- coxph(Surv(survival_time, vital_status) ~ age_at_diagnosis +
  #                    status,
  #                  data = survival.df, na.action = na.omit)
  #   }
  # }

  temp <- cox.zph(cox)

  summary(cox)

  objsurv = survfit(Surv(survival_time, vital_status) ~ status, data = survival.df)


  if (length(objsurv$strata) == 0) {
    P = ggplot()
  } else {
    P = ggsurvplot(objsurv, data = survival.df,
                   font.legend = c(14, "plain", "black"),
                   legend.title = element_blank(),
                   legend.labs = gsub("status=", "", names(objsurv$strata)),
                   # palette = "jco",
                   xlab = "Days",
                   legend = legend_pos,
                   conf.int = conf.int) +
      guides(colour = guide_legend(nrow = length(objsurv$strata)))

    P$plot = P$plot + theme(legend.background = element_rect(fill='transparent'),
                            legend.box.background =
                              element_rect(fill='transparent', linewidth = 0))
  }

  return(list(survival.df = survival.df, coxout = cox, survP = P))
}


#' Running batch survival analysis tests for a set of interactions. Returns a
#' plotlist for all the interaction test which didn't fail.
#' @param sig.sig.tissues.matrix A matrix with rows and columns with signatures,
#' and the elements are comma-separated tissue names where that interaction is
#' observed.
#' @param dataset The signature values for all samples. E.g. PCAWG.full.subset.ann
#' @param clin.df The dataframe with clinical info.
#' @param with.total.muts If TRUE the total number of mutations in the samples will
#' be provided as a confounder to the model. Default: TRUE
#' @param tmb.logged If TRUE the tumor mutational burden will be logged.
#' Default: TRUE
#' @param binary.status If TRUE, the model will compare samples with both signatures
#' with all the other samples having either of the signatures or none. Default:
#' FALSE
#' @param epistatic The survival model tests the interaction of int1 and int2 as
#' int1 * int2. Default: FALSE
#' @param legend.pos Position of the legend in the survival plot.
#' @export

get_surv_plotlist = function(sig.sig.tissues.matrix,
                             dataset,
                             clin.df,
                             with.total.muts = TRUE,
                             tmb.logged = TRUE,
                             binary.status = FALSE,
                             epistatic = FALSE,
                             legend.pos = c(0.8, 0.8)) {

  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                       base_size = 10,
                       padding = unit(c(2, 4), "mm"))

  if ( all(rownames(sig.sig.tissues.matrix) == colnames(sig.sig.tissues.matrix) ) ) {
    sig.sig.tissues.matrix[lower.tri(sig.sig.tissues.matrix, diag = TRUE)] = NA
  }
  int.indices = which(! is.na(sig.sig.tissues.matrix), arr.ind = T)


  out_plotlist = list()
  j = 0
  for (i in 1:nrow(int.indices)) {

    indices = int.indices[i, ]
    sig1 = rownames(sig.sig.tissues.matrix)[indices[1]]
    sig2 = colnames(sig.sig.tissues.matrix)[indices[2]]

    tissues = strsplit(sig.sig.tissues.matrix[indices[1], indices[2]], split = ", ")[[1]]
    for (tissue in tissues) {
      cat("Attempting cox survival for:", tissue, "::", sig1, "+", sig2, "j = ", j, "\n")
      try({surv.out = survival_for_interactions(dataset = dataset,
                                                signatures = c(sig1, sig2),
                                                tissues = tissue,
                                                clin.df = clin.df,
                                                legend_pos = legend.pos,
                                                with.total.muts = with.total.muts,
                                                tmb.logged = tmb.logged,
                                                binary.status = binary.status,
                                                epistatic = epistatic)

      cox.coefs = round(summary(surv.out$coxout)$coefficients, 2)
      rownames(cox.coefs) = gsub("status", "", rownames(cox.coefs))
      nosig = all(summary(surv.out$coxout)$coefficients[, "Pr(>|z|)"] > 0.05, na.rm = TRUE)

      if ( nosig ) {
        cat("\tNone of the covariates are significant. Skipping.\n")
        next
      }

      survp = surv.out$survP
      survp$plot = survp$plot + ggtitle(paste(tissue, "::", sig1, "+", sig2)) +
        theme(plot.title = element_text(size = 14))

      tblgrob = tableGrob(cox.coefs, theme=tt)

      coefs.ggforest = ggforest(model = surv.out$coxout, data = surv.out$survival.df)

      coefs.table = tableGrob(cox.coefs)
      # p = ggarrange(survp$plot, tblgrob, nrow = 2)
      p = cowplot::plot_grid(survp$plot, coefs.ggforest, coefs.table, nrow = 3,
                             rel_heights = c(1.5, 0.7, 0.7),
                             rel_widths = c(1, 1.5, 1.5))

      j = j + 1
      out_plotlist[[j]] = p} )
    }
  }
  return(out_plotlist)
}


#' Running batch survival analysis tests for a set of interactions. Returns a
#' list of cox regression outputs for all the interaction tests which didn't fail.
#' @param sig.sig.tissues.matrix A matrix with rows and columns with signatures,
#' and the elements are comma-separated tissue names where that interaction is
#' observed.
#' @param dataset The signature values for all samples. E.g. PCAWG.full.subset.ann
#' @param clin.df The dataframe with clinical info.
#' @param with.total.muts If TRUE the total number of mutations in the samples will
#' be provided as a confounder to the model. Default: TRUE
#' @param tmb.logged If TRUE the tumor mutational burden will be logged.
#' Default: TRUE
#' @param binary.status If TRUE, the model will compare samples with both signatures
#' with all the other samples having either of the signatures or none. Default:
#' FALSE
#' @param epistatic The survival model tests the interaction of int1 and int2 as
#' int1 * int2. Default: FALSE
#' @export

get_surv_coxlist = function(sig.sig.tissues.matrix,
                            dataset,
                            clin.df,
                            with.total.muts = TRUE,
                            tmb.logged = TRUE,
                            binary.status = FALSE,
                            epistatic = FALSE) {

  tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                       base_size = 10,
                       padding = unit(c(2, 4), "mm"))

  if ( all(rownames(sig.sig.tissues.matrix) == colnames(sig.sig.tissues.matrix) ) ) {
    sig.sig.tissues.matrix[lower.tri(sig.sig.tissues.matrix, diag = TRUE)] = NA
  }
  int.indices = which(! is.na(sig.sig.tissues.matrix), arr.ind = T)


  out_coxlist = list()
  j = 0
  for (i in 1:nrow(int.indices)) {

    indices = int.indices[i, ]
    sig1 = rownames(sig.sig.tissues.matrix)[indices[1]]
    sig2 = colnames(sig.sig.tissues.matrix)[indices[2]]

    tissues = strsplit(sig.sig.tissues.matrix[indices[1], indices[2]], split = ", ")[[1]]
    for (tissue in tissues) {
      cat("Attempting cox survival for:", tissue, "::", sig1, "+", sig2, "j = ", j, "\n")
      try({surv.out = survival_for_interactions(dataset = dataset,
                                                signatures = c(sig1, sig2),
                                                tissues = tissue,
                                                clin.df = clin.df,
                                                with.total.muts = with.total.muts,
                                                tmb.logged = tmb.logged,
                                                binary.status = binary.status)

      cox.coefs = round(summary(surv.out$coxout)$coefficients, 2)
      rownames(cox.coefs) = gsub("status", "", rownames(cox.coefs))
      nosig = all(summary(surv.out$coxout)$coefficients[, "Pr(>|z|)"] > 0.05, na.rm = TRUE)

      if ( nosig ) {
        cat("\tNone of the covariates are significant. Skipping.\n")
        next
      }


      j = j + 1
      out_coxlist[[paste(tissue, "::", sig1, "+", sig2)]] = surv.out$coxout} )
    }
  }
  return(out_coxlist)
}

#' Running batch survival analysis tests for a set of interactions. Returns a
#' list of best cox regression outputs for all the interaction tests
#' which didn't fail. The models are selected based on whether they show an
#' effect for the interaction and among those models with one with the highest
#' loglik is picked.
#' @param sig.sig.tissues.matrix A matrix with rows and columns with signatures,
#' and the elements are comma-separated tissue names where that interaction is
#' observed.
#' @param dataset The signature values for all samples. E.g. PCAWG.full.subset.ann
#' @param clin.df The dataframe with clinical info.
#' @param with.total.muts If TRUE the total number of mutations in the samples will
#' be provided as a confounder to the model. Default: TRUE
#' @param tmb.logged If TRUE the tumor mutational burden will be logged.
#' Default: TRUE
#' @param binary.status If TRUE, the model will compare samples with both signatures
#' with all the other samples having either of the signatures or none. Default:
#' FALSE

get_surv_bestcoxlist = function(sig.sig.tissues.matrix,
                                dataset,
                                clin.df,
                                with.total.muts = TRUE,
                                tmb.logged = TRUE,
                                binary.status = FALSE) {

  # tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
  #                      base_size = 10,
  #                      padding = unit(c(2, 4), "mm"))

  if ( all(rownames(sig.sig.tissues.matrix) == colnames(sig.sig.tissues.matrix) ) ) {
    sig.sig.tissues.matrix[lower.tri(sig.sig.tissues.matrix, diag = TRUE)] = NA
  }
  int.indices = which(! is.na(sig.sig.tissues.matrix), arr.ind = T)


  out_coxlist = list()
  j = 0
  for (i in 1:nrow(int.indices)) {

    indices = int.indices[i, ]
    sig1 = rownames(sig.sig.tissues.matrix)[indices[1]]
    sig2 = colnames(sig.sig.tissues.matrix)[indices[2]]

    tissues = strsplit(sig.sig.tissues.matrix[indices[1], indices[2]], split = ", ")[[1]]
    for (tissue in tissues) {
      cat("Attempting cox survival for:", tissue, "::", sig1, "+", sig2, "j = ", j, "\n")
      try({surv.out = survival_for_interactions(dataset = dataset,
                                                signatures = c(sig1, sig2),
                                                tissues = tissue,
                                                clin.df = clin.df,
                                                with.total.muts = with.total.muts,
                                                tmb.logged = tmb.logged,
                                                binary.status = binary.status)

      cox.coefs = round(summary(surv.out$coxout)$coefficients, 2)
      rownames(cox.coefs) = gsub("status", "", rownames(cox.coefs))
      nosig = all(summary(surv.out$coxout)$coefficients[, "Pr(>|z|)"] > 0.05, na.rm = TRUE)

      if ( nosig ) {
        cat("\tNone of the covariates are significant. Skipping.\n")
        next
      }


      j = j + 1
      out_coxlist[[paste(tissue, "::", sig1, "+", sig2)]] = surv.out$coxout} )
    }
  }
  return(out_coxlist)
}


#' Tests a range of models
#' @export

pick_survival_model_int = function(dataset = dataset,
                                   signatures = c(sig1, sig2),
                                   tissues = tissue,
                                   clin.df = clin.df,
                                   param.values,
                                   min.sample.fraction = 0
                                   # with.total.muts = with.total.muts,
                                   # tmb.logged = tmb.logged,
                                   # binary.status = binary.status
) {
  ### get all possible param values
  if (missing(param.values)) {
    param.values = list("with.total.muts" = c(TRUE, FALSE),
                        "tmb.logged" = c(TRUE, FALSE),
                        "binary.status" = c(TRUE, FALSE)
    )
  }
  all.param.combinations = expand.grid(param.values)

  ### exclude the impossible combinations
  # impossible.comb = c("with.total.muts" = FALSE, "tmb.logged" = TRUE)
  filtered.combinations = all.param.combinations %>%
    mutate(tmb.logged = ifelse(with.total.muts == FALSE, "FALSE", tmb.logged),
           binary.status = ifelse(epistatic, "FALSE", binary.status)) %>%
    unique()

  ### running the survival_for_interactions for given model
  lambda = function(with.total.muts, tmb.logged, binary.status, epistatic) {
    survival_for_interactions(dataset = dataset,
                              signatures = signatures,
                              tissues = tissues,
                              clin.df = clin.df,
                              with.total.muts = with.total.muts,
                              tmb.logged = tmb.logged,
                              binary.status = binary.status,
                              epistatic = epistatic)
  }

  best.model.loglik = -Inf
  best.model = NULL
  for (i in 1:nrow(filtered.combinations)) {
    param.input = filtered.combinations[i,, drop = FALSE] %>% as.list

    cat("i = ", i, "best.model.loglik = ", best.model.loglik, "\n" )

    if (param.input$epistatic) {
      param.of.interactions = paste(sort(signatures), collapse = "*")
    } else {
      param.of.interactions = paste0("status",
                                     paste(sort(signatures), collapse = "+"), collapse = "")
    }

    try({
      # print(param.input)
      ### Applying the survival_for_interactions
      test.model = do.call(lambda,  param.input)
      cat("model.loglik = ", test.model$coxout$loglik[2], "\n" )

      ### Filtering for minimum sample fraction
      sample.counts = test.model$survival.df$exists__12 %>% table()
      sample.fractions = sample.counts / sum(sample.counts)
      minority.sample.fraction = min(sample.fractions)

      if (minority.sample.fraction < min.sample.fraction) {
        cat("For Tissues == ", tissues, " and signatures == ", signatures,
            "the minimum sample threshold of ", min.sample.fraction, " has not been met.\n",
            "Skipping.\n")
        break
      }

      # return(test.model)
      model.coxout = test.model$coxout
      p.val.of.interaction = summary(test.model$coxout)$coefficients[param.of.interactions,5]
      p.val.of.interaction = ifelse(is.na(p.val.of.interaction), 1, p.val.of.interaction)

      # if (p.val.of.interaction < 0.05) {
      #   if (test.model$coxout$loglik[2] > best.model.loglik) {
      #     best.model = list(params = param.input, out.model = test.model,
      #                       minority.smp.fraction = minority.sample.fraction)
      #
      #     best.model.loglik = test.model$loglik
      #     # cat("\ntissue = ", tissues, "\n")
      #     # print(unlist(param.input) )
      #     # cat("\n")
      #   }
      # }

      if(test.model$coxout$loglik[2] > best.model.loglik) {
        if ( is.null(best.model) ) {
          best.model = list(params = param.input, out.model = test.model,
                            minority.smp.fraction = minority.sample.fraction, ind = i)
          best.model.loglik = test.model$coxout$loglik[2]
        } else {
          cat("PLR test for models", i, "and ", best.model$ind, "\n")
          plrtest.out = plrtest(test.model$coxout, best.model$out.model$coxout,
                                nested = FALSE)

          if (plrtest.out$pLRTAB < 0.05) {
            cat ("Model1 and model2 are distinguisable!!!!######!!!!!")
            if (plrtest.out$pLRTA < 0.05) {
              cat("Model1 fits better according to PLR.\n")
              best.model = list(params = param.input, out.model = test.model,
                                minority.smp.fraction = minority.sample.fraction, ind = i)
              best.model.loglik = test.model$loglik[2]
            } else {
              if (test.model$coxout$loglik[2] > best.model.loglik ) {
                best.model = list(params = param.input, out.model = test.model,
                                  minority.smp.fraction = minority.sample.fraction, ind = i)
                best.model.loglik = test.model$loglik[2]
              }
            }
          } else {
            cat ("Model1 and model2 are not distinguisable!!!!.....!!!!!")
          }
        }
        }
    } )
  }
  return(best.model)
}


#' Running batch survival analysis tests for a set of interactions. Returns a
#' list of cox regression outputs with the best selected model. If the tissue
#' supports a survival effect for a given interaction, then this model is selected.
#' Among several models, the one with the lowest loglik is .
#' @param sig.sig.tissues.matrix A matrix with rows and columns with signatures,
#' and the elements are comma-separated tissue names where that interaction is
#' observed.
#' @param dataset The signature values for all samples. E.g. PCAWG.full.subset.ann
#' @param clin.df The dataframe with clinical info.
#' @param with.total.muts If TRUE the total number of mutations in the samples will
#' be provided as a confounder to the model. Default: TRUE
#' @param tmb.logged If TRUE the tumor mutational burden will be logged.
#' Default: TRUE
#' @param binary.status If TRUE, the model will compare samples with both signatures
#' with all the other samples having either of the signatures or none. Default:
#' FALSE
#' @export

get_surv_best_model = function(sig.sig.tissues.matrix,
                               dataset,
                               clin.df,
                               param.list,
                               min.sample.fraction = 0) {

  tt <- gridExtra::ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                                  base_size = 10,
                                  padding = unit(c(2, 4), "mm"))

  if (missing(param.list)) {
    param.list = list("with.total.muts" = c(TRUE, FALSE),
                      "tmb.logged" = c(TRUE, FALSE),
                      "binary.status" = c(TRUE, FALSE))
  }

  if ( all(rownames(sig.sig.tissues.matrix) == colnames(sig.sig.tissues.matrix) ) ) {
    sig.sig.tissues.matrix[lower.tri(sig.sig.tissues.matrix, diag = TRUE)] = NA
  }
  int.indices = which(! is.na(sig.sig.tissues.matrix), arr.ind = T)

  out_coxlist = list()
  j = 0
  for (i in 1:nrow(int.indices)) {

    indices = int.indices[i, ]
    sig1 = rownames(sig.sig.tissues.matrix)[indices[1]]
    sig2 = colnames(sig.sig.tissues.matrix)[indices[2]]

    tissues = strsplit(sig.sig.tissues.matrix[indices[1], indices[2]], split = ", ")[[1]]
    for (tissue in tissues) {
      cat("Attempting cox survival for:", tissue, "::", sig1, "+", sig2, "j = ", j, "\n")
      try({surv.out = pick_survival_model_int(dataset = dataset,
                                              signatures = c(sig1, sig2),
                                              tissues = tissue,
                                              clin.df = clin.df,
                                              param.values = param.list,
                                              min.sample.fraction = min.sample.fraction)
      if ( is.null(surv.out$out.model) ) {
        cat("\tThe interaction is not significant. Skipping.\n")
        next
      }

      j = j + 1
      out_coxlist[[paste(tissue, "::", sig1, "+", sig2)]] =
        list(input.params = unlist(surv.out$params),
             model = surv.out$out.model$coxout,
             minority.smp.fraction = surv.out$minority.smp.fraction)} )
    }
  }
  return(out_coxlist)
}





#' A function to generate default colors in ggplots.
#' @param n number of colors.
#'
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' The function outputs mutational profiles for mutations in genes
#' in the pathway of interest and for all other mutations.
#' @param pathway One of the 10 oncogenic pathway names.
#' @param tissue The tissues where the profiles have to be assessed.
#' @param pathways.df The data.frame with sample names and pathway alterations.
#' @param vcf.dir Directory to vcf files
#' @return A list with mutational profiles for all mutations and for those in
#' the pathway of interest.
#' @export

pcawg_pathway_profiles = function(pathway, tissue, pathways.df, vcf.dir =
                                    here("data/raw/PCAWG/vcfs/consensus_snv_indel/snv_mnv/")) {
  path.genes = all.gene.path %>% filter(Pathway == pathway) %>% pull(Gene)

  path.gene.symbols <- biomaRt::select(org.Hs.eg.db,
                                       keys = path.genes,
                                       columns = c("SYMBOL","ENTREZID"),
                                       keytype = "SYMBOL")

  path.gene.coords <- biomaRt::select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      keys = path.gene.symbols$ENTREZID,
                                      columns = c("GENEID", "TXID", "TXCHROM", "TXSTART", "TXEND"),
                                      keytype = "GENEID")

  path_gene_granges = GRanges(seqnames = path.gene.coords$TXCHROM,
                              IRanges(start = path.gene.coords$TXSTART,
                                      end = path.gene.coords$TXEND)) %>%
    disjoin()

  tissue.samples = pathways.df %>% filter(Cancer.Types == tissue) %>%
    pull(sample_id)

  vcf.files <- list.files(vcf.dir,
                          pattern = paste(tissue.samples, collapse = "|"), full.names = TRUE
  )
  vcf.files = vcf.files[!grepl("tbi", vcf.files)]

  vcf.file.order = gsub(".consensus.*$", "", basename(vcf.files))
  vcf.samples = pathways.df[match(vcf.file.order, pathways.df$sample_id),] %>%
    pull(donor_id)

  cat("Reading vcf files... This may take a while.\n")
  grl <- MutationalPatterns::read_vcfs_as_granges(vcf.files, vcf.samples, "BSgenome.Hsapiens.UCSC.hg19")
  cat("Done.\n")

  grl.path = lapply(grl, function(x)
    subsetByOverlaps(x, path_gene_granges) )

  type.occurrences.path <- mut_type_occurrences(grl.path, ref.genome)

  # total.type.occurrences = colSums(type.occurrences.path)


  type.occurrences <- mut_type_occurrences(grl, ref.genome)
  return(list(path.profiles = type.occurrences.path,
              all.profiles = type.occurrences))
}

#' Hypergeometric test for C>T mutations at CpG sites between two profiles.
#' @param seven_channel_profile1 First profile
#' @param seven_channel_profile2 Second profile
#' @param lower.tail Parameter will be passed to phyper
#' @return P value or 1-P of a hypergeometric test.
#' @export

test_CpG = function(seven_channel_profile1, seven_channel_profile2, lower.tail = TRUE) {

  ### Following the logic of the answer
  ### https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper

  p1 = seven_channel_profile1
  p2 = seven_channel_profile2

  q = p1["C>T at CpG"]
  m = p1["C>T at CpG"] + p2["C>T at CpG"]
  n = p1["C>T"] + p2["C>T"] - m
  k = p1["C>T"]

  return(phyper(q, m, n, k, lower.tail))
}



#' Hypergeometric test for C>T mutations at CpG sites between two profiles.
#' @param seven_channel_profile1 First profile
#' @param seven_channel_profile2 Second profile
#' @param lower.tail Parameter will be passed to phyper
#' @return P value or 1-P of a hypergeometric test.

test_mutation = function(mutation, seven_channel_profile1, seven_channel_profile2, lower.tail = TRUE) {

  ### Following the logic of the answer
  ### https://stackoverflow.com/questions/8382806/hypergeometric-test-phyper

  p1 = seven_channel_profile1
  p2 = seven_channel_profile2

  total.p1 = p1["C>A"] + p1["C>G"] + p1["C>T"] + p1["T>A"] + p1["T>C"] + p1["T>G"]
  total.p2 = p2["C>A"] + p2["C>G"] + p2["C>T"] + p2["T>A"] + p2["T>C"] + p2["T>G"]


  q = p1[mutation]
  m = p1[mutation] + p2[mutation]
  n = p1[total.p1] + p2[total.p2] - m
  k = total.p1

  return(phyper(q, m, n, k, lower.tail))
}


#' A wrapper to call GES scoring tool from within R.
#'
#' @param mat Input matrix or a data frame - the feature matrix
#' @param tmpdir Directory where the temporary files should be written.
#'
#' @details Initial implementation in git repo:
#' Biwei-Huang/Generalized-Score-Functions-for-Causal-Discovery/
#' G(i,j) = 1: i->j (the edge is from i to j)
#' G(i,j) = -1: i-j (the direction between i and j is not determined)
#' G(i,j) = 0:  i j (no causal edge between i and j)
#'
#' @return DAG adjancency matrix from GES method.
#'

run_GES_octave_wrapper = function(mat, tmpdir = here("tmp")) {

  tmp.filename = stringi::stri_rand_strings(1, 10, '[a-zA-Z]')
  tmp.in = paste0(file.path(tmpdir, tmp.filename), ".in")
  tmp.out = paste0(file.path(tmpdir, tmp.filename), ".out")

  if (is.matrix(mat)) {
    mat = as.data.frame(mat)
  }

  write.table(mat, file = tmp.in, sep = "\t",
              quote = FALSE, col.names = TRUE, row.names = FALSE)
  cat("Running octave...\n")
  octave_out = system(paste("octave ", here("octave/run_melanoma.m"),
                            tmp.in, tmp.out) )

  if (octave_out) {
    stop("Octave encountered an error.")
  }

  GES = read.table(tmp.out, sep = "\t",
                   col.names = colnames(mat),
                   row.names = colnames(mat))
  return(GES)
}



#' Order matrix/dataframe columns/rows in increasing order
#' @param input.mat input data
#' @export

order_matrix_rc = function(input.mat) {

  cnames = sort(colnames(input.mat))
  rnames = sort(rownames(input.mat))

  return(input.mat[rnames, cnames])

}


#' Plots a signature heatmap with pathway annotations
#'
#' @param tissue Tissue
#' @param signatures Signatures dataframe, signatures start from column 4.
#' Column 2 is Sample.Names.
#' @param pathways Pathway mutations dataframe, pathways start from column 4,
#' contains donor_id.
#' @param border_color Border color of cells. This parameter is there to control
#' removing border color with NA. Pheatmap doesn't properly remove it
#' and the border is still present when saving the plot. Default: gray60.
#' Should be NA to remove the border_color.
#' @param ... Params passed to pheatmap.
#' @export

pathways_signatures_heatmap = function(tissue, signatures, pathways,
                                       border_color = "grey60", ...) {

  tissue.sigs = subset_tissue(signatures, tissue)

  common.samples = intersect(pathways$donor_id, tissue.sigs$Sample.Names)

  tissue.pathways = pathways %>% filter(donor_id %in% common.samples) %>% arrange(match(donor_id, common.samples))

  tissue.pathways = cbind(tissue.pathways["donor_id"],
                          tissue.pathways[, names(which(colSums(
                            tissue.pathways[4:ncol(tissue.pathways)], na.rm = TRUE) > 0))])

  tissue.sigs =  tissue.sigs %>% filter(Sample.Names %in% common.samples) %>%
    arrange(match(Sample.Names, common.samples))

  rownames(tissue.sigs) = tissue.sigs$Sample.Names

  pp = pheatmap(log(tissue.sigs[4:ncol(tissue.sigs)] + 1), annotation_row = tissue.pathways[, 2:ncol(tissue.pathways)],
                annotation_legend = FALSE, color = viridis(15),
                show_rownames = FALSE, ...)


  if (is.na(border_color)) {
    # https://stackoverflow.com/questions/44318690/no-border-color-in-pheatmap-in-r
    grob_classes <- purrr::map(pp$gtable$grobs, class)

    ### Finding the histogram

    idx_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'gTree' %in% cl))[1]
    grob_names <- names(pp$gtable$grobs[[idx_grob]]$children)
    idx_rect <- grob_names[grep('rect', grob_names)][1]

    ## Remove borders around cells
    pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$col <- pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$fill
    pp$gtable$grobs[[idx_grob]]$children[[idx_rect]]$gp$lwd <- 0

    ### Finding annotation columns to the left

    annot_grob <- which(purrr::map_lgl(grob_classes, function(cl) 'rect' %in% cl))[1]

    ## Remove borders around cells
    pp$gtable$grobs[[annot_grob]]$gp$col = NA
  }
  return(pp)
}


#' Summarizes coxph statistics from a list of coxph objects. List names are the
#' names of interactions inside the tissues.
#' @param interaction.summaries the list of coxph objects
#' @param type if not NULL, the output has one more column with this variable.
#' Can be used to indicate if the interaction is positive or negative.
#' @return A dataframe summarizing the regression output for all the
#' models and parameters.
#' @export

HR_summary_for_all = function(interaction.summaries, type = NULL) {

  get_surv_vector = function(cox.out, cond, model.params = NULL) {
    summary.out = summary(cox.out)
    coeffs = summary.out$coefficients
    conf.int = summary.out$conf.int

    P.val = coeffs[, 5]

    sig.star = sapply(P.val, get_sig_stars)

    lower.95 = conf.int[, 3]
    upper.95 = conf.int[, 4]
    estimate = conf.int[, 1]
    out.data = data.frame(params = gsub(pattern = "status", "", rownames(conf.int) ),
                          cond = cond,
                          estimate = estimate,
                          lower.95 = lower.95,
                          upper.95 = upper.95,
                          P.val = P.val,
                          sig.star = sig.star)

    if (!is.null(model.params)) {
      out.data = cbind(out.data, model.params %>% t() %>%
                         as.data.frame() %>%
                         slice(rep(1:n(), each = nrow(out.data) ) ) )
    }
    return(out.data)
  }

  if (all(c("input.params", "model") %in% names(interaction.summaries[[1]])) ) {
    models.list = lapply(interaction.summaries, function(x) x[["model"]])
    model.params.list = lapply(interaction.summaries, function(x) x[["input.params"]])
  } else {
    models.list = interaction.summaries
    model.params.list = NULL
  }

  all.cond.summaries.list = lapply(names(models.list), function(x)
    get_surv_vector(models.list[[x]], x, model.params.list[[x]]))
  all.params.ever = do.call(rbind, all.cond.summaries.list)
  all.params.ever$tissue = sapply(strsplit(all.params.ever$cond, " :: "), function(x) x[1])
  all.params.ever$int = sapply(strsplit(all.params.ever$cond, " :: "), function(x) x[2])

  if (! is.null(type)) {
    all.params.ever$type = type
  }

  return(all.params.ever)
}

#' (plot hazard ratio for variables)
#' Makes a forest plot for hazard ratios for given parameter
#' @param all.conds.df Data.frame with all the parameters, hazard ratios,
#' confidence intervals. Generated with function HR_summary_for_all.
#' @param param Parameter to be plotted
#' @param average If TRUE the parameter values for a given tissues should be
#' averaged. Default:TRUE.
#' @param log.HR If TRUE the Hazard ratio will be logged. Default: FALSE.
#' @export

plot_HR_vars = function(all.conds.df, param, average = TRUE, log.HR = FALSE, no_stripes = TRUE) {

  tissue.increasing.order = all.conds.df %>%
    filter(params == param) %>%
    arrange(estimate) %>%
    pull(tissue) %>%
    unique()

  all.conds.df = all.conds.df %>%
    filter(params == param) %>%
    arrange(estimate) %>%
    mutate(tissue = factor(tissue, levels = tissue.increasing.order))
  # all.conds.df$tissue = factor(all.conds.df$tissue,
  #                       levels = sort(unique(all.conds.df$tissue)))
  if (average) {
    df = all.conds.df %>%
      group_by(tissue) %>%
      summarize(estimate = mean(estimate),
                lower.95 = mean(lower.95),
                upper.95 = mean(upper.95),
                P.val = mean(P.val),
                sig.star = get_sig_stars(P.val)) %>%
      mutate(sig.effect = as.character(sign((upper.95 - 1) * (lower.95 - 1) ) ))
    if (log.HR) {
      df = df %>% mutate(estimate = log(estimate),
                         lower.95 = log(lower.95),
                         upper.95 = log(upper.95),
                         sig.effect = as.character(sign(upper.95 * lower.95) ))
    }

    p = ggplot(df, aes(y = tissue, x = estimate, alpha = sig.effect) )

  } else {
    df = all.conds.df %>%
      mutate(sig.effect = as.character(sign((upper.95 - 1) * (lower.95 - 1) ) ) )
    if (log.HR) {
      df = df %>% mutate(estimate = log(estimate),
                         lower.95 = log(lower.95),
                         upper.95 = log(upper.95),
                         sig.effect = as.character(sign(upper.95 * lower.95) ) )
      df = df %>% mutate(sig.effect = as.character(sign(upper.95 * lower.95) ) )
      p = ggplot(df, aes(y = tissue, x = estimate, color = int, alpha = sig.effect) )
    }
  }

  p = p + scale_alpha_manual(values = c("-1" = 0.4, "1" = 1))


  p = p + geom_point(position=position_dodge(1), shape = 15, size = 3) +
    geom_errorbar(aes(y = tissue, xmin = lower.95, xmax = upper.95),
                  width=.2,
                  position=position_dodge(1)) +
    scale_y_discrete(position = "right")

  if (log.HR) {
    p = p + geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
      xlab("log(Hazard Ratio)")
  } else {
    p = p + geom_vline(xintercept = 1, linetype="dashed", color = "gray") +
      xlab("Hazard Ratio")
  }

  tissue.mapping = all.conds.df%>% select(tissue, data) %>% unique()
  tissue.map = setNames(substr(tissue.mapping$data, 1, 1), tissue.mapping$tissue)

  pgbuild = ggplot_build(p)
  xlimits = pgbuild$layout$panel_scales_x[[1]]$range$range
  ylimits = pgbuild$layout$panel_scales_y[[1]]$range$range

  annot.shift = (xlimits[2] - xlimits[1] )/ 10

  p = p +
    geom_hline(yintercept = (seq(ylimits) - 0.5 )[2:length(ylimits) ],
               color = "gray90")  +
    annotate("text", x = xlimits[1] - annot.shift, y = ylimits, label = tissue.map[ylimits]) +
    coord_cartesian(xlim = xlimits, clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 0, 0.7, "cm"))

  if (log.HR) {
    p = p + xlab("log(HR)")
  } else {
    p = p + xlab("HR")
  }

  return(p)
}


#' Creates summaries from survival summaries generated by HR_summary_for_all
#' @param data summary output created by HR_summary_for_all
#' @param param the name of the parameter to summarize
#' @export

plot_param_piechart = function(data, param) {
  data.tissues = data %>%
    filter(params == param) %>%
    mutate(HR.sign = estimate > 1,
           significance = ifelse(sig.star != " ", TRUE, FALSE),
           HR.sign = ifelse(significance == FALSE, FALSE, HR.sign)) %>%
    group_by(tissue, significance, HR.sign) %>%
    select(tissue, significance, HR.sign) %>% unique() %>%
    arrange(-significance, -HR.sign)

  data.tissues = data.tissues[ !(duplicated(data.tissues$tissue)),]

  p = data.tissues %>%
    group_by(significance, HR.sign) %>%
    summarize(counts = n()) %>%
    mutate(colorcode = ifelse(significance == FALSE, "no effect",
                              ifelse(HR.sign, "positive", "negative"))) %>%
    ggplot(aes(x = "", y = counts, fill = colorcode)) +
    geom_col(color = "white") +
    geom_text(aes(label = counts),
              position = position_stack(vjust = 0.5), size =7) +
    scale_fill_manual(values = c(`no effect` = "gray90",positive = "coral1", negative = "deepskyblue3"),
                      name = "Effect on HR") +
    coord_polar(theta = "y") + theme_void()
  # theme(legend.position = "none")
  return(p)
}


#' Make forest plot for all the significant interaction paramters from the
#' survival models.
#' @param data The summary data across tissues and
#' models generated with HR_summary_for_all.
#' @param log.HR Indicated whether the hazard ratio is logged or not.
#' @export

plot_sigint_forest = function(data, log.HR = TRUE) {
  data.processed = data %>%
    drop_na() %>%
    filter(P.val < 0.05, ! is.infinite(lower.95), ! is.infinite(upper.95),
           ! is.infinite(log(lower.95) ), ! is.infinite(log(upper.95)) ) %>%
    filter(grepl("+", params,fixed = TRUE)) %>%
    filter(! params %in% c("age_at_diagnosis", "log(total_muts + 1)"))

  param.order = data.processed %>%
    arrange(tissue, estimate) %>%
    pull(params) %>% unique() %>% rev()

  data.processed = data.processed %>%
    mutate(params = factor(params, levels = param.order))

  pp = data.processed %>%
    ggplot(aes(y = params, x = log(estimate), color = tissue, shape = type)) +
    geom_point(position=position_dodge(0.7), size = 3) +
    geom_errorbar(aes(y = params, xmin = log(lower.95), xmax = log(upper.95)),
                  width=.2,
                  position=position_dodge(0.7)) +
    scale_y_discrete(position = "right")

  pgbuild = ggplot_build(pp)
  xlimits = pgbuild$layout$panel_scales_x[[1]]$range$range
  ylimits = pgbuild$layout$panel_scales_y[[1]]$range$range

  pp = pp +
    geom_hline(yintercept = (seq(ylimits) - 0.5 )[2:length(ylimits) ], alpha = 0.1) +
    scale_shape_manual(values = c("pos" = 17, "neg" = 19),
                       labels = c("pos" = "positive", "neg" = "negative"),
                       name = "Interaction type") +
    scale_color_discrete(name = "Tissue")

  if (log.HR) {
    pp = pp +
      geom_vline(xintercept = 0, linetype="dashed", color = "gray") +
      xlab("log(HR)")
  } else {
    pp = pp +
      geom_vline(xintercept = 1, linetype="dashed", color = "gray") +
      xlab("HR")
  }

  pp = pp + ylab("") +
    # coord_cartesian(xlim = c(-1, 20)) +
    theme_classic(base_size = 15) +
    theme(axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  return(pp)
}

#' Plot interaction network with a mixed layout, where the clock-like
#' signatures are in the middle and all the other signatures are on a circular
#' layout.
#' @param graph.input The network of interactions as a tbl_graph
#' @param central.nodes The nodes that should be in the center.
#' Default: c("Ageing" , "SBS5", "SBS40")
#' @param circular.node.order Can be used to specify the order of the nodes in the
#' outer layer. Default: NULL
#' @param edge.width.breaks Can be used to manually specify the edge widths.
#' Default: NULL
#' @param lbm legend.box.margin to be passed to theme
#' @details Beware!!! The node colors are predefined for different etiology groups.
#' The etiology groups are derived from signature.annotations object
#' @return a ggraph plot of the interaction network.
#' @export

plot_mixed_layout = function(graph.input,
                             central.nodes = c("Ageing" , "SBS5", "SBS40"),
                             circular.node.order = NULL, edge.width.breaks = NULL,
                             lbm = margin(0,0,0, 1.1, unit = "in")) {


  all.nodes = graph.input %>%
    activate(nodes) %>%
    as.data.frame() %>%
    pull(name)

  graph.input = graph.input %>%
    activate(nodes) %>%
    mutate(annot.class = signature.annotations %>%
             filter(Annotation %in% all.nodes) %>%
             select(Origin, Annotation) %>%
             unique() %>%
             arrange(factor(Annotation, levels = all.nodes)) %>%
             pull(Origin))

  circle.nodes = setdiff(all.nodes, central.nodes)
  #
  # lay1 = layout_in_circle(induced_subgraph(graph.input, circle.nodes))
  # lay2 = layout_nicely(induced_subgraph(graph.input, central.nodes))
  # lay3 <- rbind(lay1 + 2, lay2 *  10)
  #

  central.subgraph = graph.input %>%
    to_subgraph(name %in% central.nodes, subset_by = "nodes")

  circle.subgraph = graph.input %>%
    to_subgraph(name %in% circle.nodes, subset_by = "nodes")

  lay.central = create_layout(central.subgraph$subgraph,
                              layout = "stress")

  if (is.null(circular.node.order)) {
    lay.circle = create_layout(circle.subgraph$subgraph,
                               layout = "circle")
  } else {
    lay.circle = create_layout(circle.subgraph$subgraph,
                               layout = "circle", order = circular.node.order)
  }

  # plot.igraph(pcawg.sig.sig.net, layout=lay3, vertex.size=8)

  node.coordinates = rbind(lay.central, lay.circle) %>% as.data.frame() %>%
    arrange(factor(name, levels = all.nodes))
  pp = graph.input %>% ggraph(x = node.coordinates$x, y = node.coordinates$y) +
    geom_edge_parallel (aes(color = int.type, width = count), alpha = 0.6) +
    scale_edge_width_continuous(breaks = edge.width.breaks, range = c(0.3,3)) +
    geom_node_point(aes(fill = annot.class),
                    color = "gray30",
                    size = 4,
                    shape = 21,
                    stroke = 0.6) +
    geom_node_label(aes(label = name#, fill = type
    ),
    color = "gray10",
    label.size = 0,
    hjust = ifelse(node.coordinates[,1] > 0, -0.25, 1.25),
    vjust = ifelse(node.coordinates[,2] > 0, -0.3, 1.15),
    # hjust = node.coordinates[ ,1] * 1.1,
    # vjust = node.coordinates[ ,2] * 1.1,
    label.r = unit(0.2, "lines"),
    label.padding = unit(0, "lines"),
    nudge_y = 0) +
    scale_edge_color_manual(values = c("positive" = "red4",
                                       "negative" = "dodgerblue3")) +
    scale_fill_manual(values = c("Endogenous" = "#CAE7B9",
                                 "Environmental" = "#f3de8a",
                                 # "Clock-like" = "#ff934f",
                                 "Clock-like" = "white",
                                 "Unknown" = "#4e6151") ) +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=lbm)
  return(pp)
}

#' Removes identical edges between two nodes
#' @param input.graph Input graph
#' @param summary.col This variable indicates the column, which should be summarized.
#' When specified, all the edges between the same nodes, which have different values
#' for this column are summarized, a new edge variable `count` is added and
#' the values of the specified column are collapsed and stored in the column called
#' `summedList`
#' @export

graph_unique_edges = function(input.graph, summary.col = NULL) {
  summary.col = "tissue"

  new.edges = input.graph %>%
    activate(edges) %>%
    as.data.frame()

  new.edges$minmax = purrr::pmap(new.edges, ~ paste0(sort(c(...)), collapse = "-") )  %>%
    unlist()

  new.edges = new.edges %>% distinct(minmax, .keep_all = TRUE) %>%
    select(-minmax) %>%
    as_tibble()

  if (!is.null(summary.col)) {
    new.edges = new.edges %>%
      group_by(across(c(-!!sym(summary.col)))) %>%
      dplyr::summarize(count = n(),
                       summedList = toString(sort(unique(!!sym(summary.col) ) ) ) ,
      ) %>%
      ungroup()
  }


  output.graph = tbl_graph(nodes = input.graph %>% activate(nodes) %>% as.data.frame(),
                           edges = new.edges,
                           directed = FALSE)

  return(output.graph)
}

#' interactive network plot
#' @param pp ggraph output from plot_mixed_layout. In fact any layout would work.
#' @param width.svg will be passed to girafe function as width_svg
#' @param height.svg will be passed to girafe function as height_svg
#' @export

plot_layout_interactive = function(pp, width.svg = NULL, height.svg = NULL) {

  ggbuild = ggplot_build(pp)

  edge.connections = get_edges()(pp$data)

  edge.connections.agg = edge.connections %>%
    group_by(x, y, xend, yend) %>%
    mutate(edge.label = paste0( paste0(int.type, ":", summedList),
                                collapse = "\n"))

  pp.int = pp +
    geom_point_interactive(aes(x, y),
                           data = pp$data, size = 5, alpha = 0) +
    geom_segment_interactive(aes(x = x, y = y, xend = xend, yend = yend,
                                 tooltip = edge.label),
                             size = 2, alpha = 0.1, color = "gray60",
                             data = edge.connections.agg)

  # pp.int = ggplot() +
  #   geom_point_interactive(aes(x = x, y = y), data = pp$data, size = 5,
  #                          tooltip = TRUE) +
  #  interactive_segments_grob(aes(x0 = x, y0 = y, x1 = xend, y1 = yend, tooltip = "blah"),
  #                            data = edge.connections)

  ppiraph = girafe(ggobj = pp.int,
                   options = list(
                     # opts_hover(css = girafe_css(
                     #     css = "fill:orange; stroke:red",
                     #     text = "stroke:blue; font-size: larger",
                     #     line = "fill:black; stroke-width:3px; opacity:1",
                     #     area = "stroke-width:3px",
                     #     point = "stroke-width:3px" ) ),
                     opts_tooltip("background-color:gray;color:white;
                                 font-style:italic;padding:6px;
                                 border-radius:7px;")
                   ),
                   width_svg = width.svg,
                   height_svg = height.svg)

  # girafe(ggobj = pp.int,
  #                options = list(
  #                      # opts_hover(css = "fill:green; stroke: blue;"),
  #                        # opts_hover(css = "fill:wheat;stroke:orange;r:5pt;"),
  #                        opts_hover(css = girafe_css(
  #                              css = "fill:orange;stroke:gray;",
  #                              text = "stroke:blue; font-size: larger",
  #                              line = "fill:none; stroke:blue",
  #                              area = "stroke-width:3px;stroke:blue",
  #                              point = "stroke-width:3px",
  #                              image = "outline:2px red"
  #                          ) ),
  #                        opts_tooltip("background-color:gray;color:white;
  #                                font-style:bold;padding:6px;
  #                                border-radius:7px;")
  #                  ))

  return(ppiraph)
}

#' Get the percentile of a value v in a vector x.
#' @param x The vector
#' @param v The value
#'
#' @export

percentile <- function(x, v) {
  sum(x <= v) / length(x)
}
