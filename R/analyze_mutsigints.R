# ### Getting all the functions used in mutsigints to compile them in this package
#
# ### Using the code from:
# # https://community.rstudio.com/t/is-there-a-way-to-extract-the-names-of-all-functions-in-an-r-script/51905/3
#
# get_calls <- function(filepath) {
#   code <- parse(filepath)
#   tokens <- as.list(code)
#   calls <- c()
#   while (TRUE) {
#     any_unpacked <- FALSE
#     for (ii in seq_along(tokens)) {
#       part <- tokens[[ii]]
#       # Calls always have the function name as the first element
#       if (missing(part)) {
#         next
#       }
#       if (is.call(part)) {
#         fun_token <- part[[1]]
#         calls <- c(calls, deparse(fun_token))
#       }
#       # Expressions have a length
#       if (length(part) > 1) {
#         tokens[[ii]] <- as.list(part)
#         any_unpacked <- TRUE
#       }
#     }
#     tokens <- unlist(tokens)
#     if (!any_unpacked) break
#   }
#   unique(calls)
# }
#
# mutsigints.analysis.Rmd = list.files("../mutsigints/analysis/", pattern = "*.Rmd", full.names = TRUE)
#
# mutsigints.analysis.R = list.files("../mutsigints/analysis/", pattern = "*.R$", full.names = TRUE)
#
# mutsigints.R.R = list.files("../mutsigints/R/", pattern = "*.R$", full.names = TRUE)
#
# mutsigints.R.R = setdiff(mutsigints.R.R, "../mutsigints/R//functions.R")
#
# all.calls = c()
#
# for (file in mutsigints.analysis.Rmd) {
#
#   output.file = file.path("../tmp", gsub(".Rmd", ".R", basename(file)))
#   knitr::purl(input = file, output = output.file)
#
#   file.calls = get_calls(output.file)
#   all.calls = unique(c(all.calls, file.calls))
# }
#
# for(file in c(mutsigints.analysis.R, mutsigints.R.R)) {
#   all.calls = unique(c(all.calls, get_calls(file)))
# }
#
# exclude.package.list = c("base", "ggplot2","stats", "dplyr", "utils", "MASS","tidygraph", "igraph", "tidyr", "ggraph" , "survminer",
#                          "devtools", "RTCGA", "TCGAbiolinks", "here")
#
#
# for(package in exclude.package.list) {
#   all.calls = setdiff(all.calls, ls(paste0("package:", package)))
# }
#
#
