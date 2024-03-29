% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{sig_pathway_int}
\alias{sig_pathway_int}
\title{Assess signature-pathway interactions across tissues for a custom function}
\usage{
sig_pathway_int(
  sigs.input,
  pathways.input,
  interaction_function,
  path.min.tissues = 30,
  p.val.threshold = 0.1,
  p.adjust = TRUE,
  method = "BH",
  ...
)
}
\arguments{
\item{sigs.input}{Signature activities. The rows correspond to samples, the
first three columns are Cancer.Types, Sample.Names, Accuracy, all the other
columns correspond to signature names. The values are signature activities.}

\item{pathways.input}{Pathway status - formatted like mutated.pathways.tissues
Pathway activities start at column 4. The first three columns are sample_id,
donor_id, Cancer.Types. Other columns correspond to pathway names. The values
correspond to number of mutations in the pathway.}

\item{interaction_function}{The function defining the metric. E.g. get_sig_path_assocs}

\item{path.min.tissues}{Minimal number of samples in each tissue to be considered}

\item{p.val.threshold}{p-value threshold for BH correction. Default: 0.05}

\item{p.adjust}{Controls if p-values should be adjusted. Default: TRUE}

\item{method}{P-value adjustement methods. Default: BH}
}
\value{
A list with the length of abundant tissues in the datasets,
where each element is the interaction matrix
}
\description{
Assess signature-pathway interactions across tissues for a custom function
}
