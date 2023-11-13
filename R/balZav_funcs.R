#' Getting symmetric pivot coordinates from a strictly positive matrix
#' Used to calculate cor_coda with rand.add = TRUE

balZav <- function(x) {
    D <- ncol(x)
    Z.av <- matrix(NA, ncol = 2, nrow = nrow(x))
    p1 <- sqrt(D - 1 + sqrt(D * (D - 2)))/sqrt(2 * D)
    if (D == 3) {
        p2 <- x[, 3]
    }
    else {
        p2 <- apply(x[, 3:D], 1, prod)
    }
    p3 <- (sqrt(D - 2) + sqrt(D))/(sqrt(D - 2) * (D - 1 + 
                                                  sqrt(D * (D - 2))))
    p4 <- 1/(D - 1 + sqrt(D * (D - 2)))
    Z.av[, 1] <- p1 * (log(x[, 1]/(x[, 2]^p4 * p2^p3)))
    Z.av[, 2] <- p1 * (log(x[, 2]/(x[, 1]^p4 * p2^p3)))

    return(Z.av)
}

#' Getting symmetric pivot coordinates from a matrix containing 0-values
#' symmetric coordinates are calculated for one sample at a time using only
#' strictly positive signatures.
#' Used to calculate cor_coda with rand.add = FALSE


balZavRow <- function(x) {
    D <- length(x)
    Z.av <- c(0,0)
    p1 <- sqrt(D - 1 + sqrt(D * (D - 2)))/sqrt(2 * D)
    if (D == 3) {
        p2 <- x[3]
    }
    else {
        p2 <- prod(x[3:D])
    }
    p3 <- (sqrt(D - 2) + sqrt(D))/(sqrt(D - 2) * (D - 1 + 
                                                  sqrt(D * (D - 2))))
    p4 <- 1/(D - 1 + sqrt(D * (D - 2)))
    Z.av[1] <- p1 * (log(x[1]/(x[2]^p4 * p2^p3)))
    Z.av[2] <- p1 * (log(x[2]/(x[1]^p4 * p2^p3)))
    return(Z.av)
}
