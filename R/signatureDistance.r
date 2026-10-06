#' Signature Distance
#'
#' This function computes the similarity between columns of a data matrix
#' @param dset1 Dataset of any type in matrix format, with features in rows and samples in columns
#' @param dset2 Optional Dataset. If provided, distance between columns of dset and dset2 are computed and reported as rows and columns, respectively; if not, distance between all possible pairs of columns from dset are computed
#' @param nn Optional size for the signature, default is either the full signature or 10 percent of it, depending or whether \code{ws}=0 or not
#' @param groups Optional vector indicating the group ID of the samples
#' @param scale. Logical, whether the data should be scaled
#' @param two.tails Logical, whether a two tails, instead of 1 tail test should be performed
#' @param ws Number indicating the exponent for the weighting the signatures, the default of 0 is uniform weighting, 1 is weighting by SD
#' @return Object of class \code{signatureDistance} as a matrix of normalized enrichment scores
#' @examples
#' data(bcellViper, package="bcellViper")
#' dd <- signatureDistance(exprs(dset))
#' dd[1:5, 1:5]
#' scale(dd)[1:5, 1:5]
#' as.matrix(as.dist(dd))[1:5, 1:5]
#' @export

signatureDistance <- function(dset1, dset2=NULL, nn=NULL, groups=NULL, scale.=TRUE, two.tails=TRUE, ws=2) {
    if (is.null(nn)) {
        if (ws>0) nn <- nrow(dset1)
        else nn <- nrow(dset1)/10
    }
    nn <- round(nn/2)
    d1 <- dset1
    if (is.null(dset2)) {
        if (ws==0 & nrow(d1)<(nn*3)) stop("Error, insufficient number of features to perform analysis", call.=FALSE)
        if (!is.null(groups)) d1 <- scaleGroups(dset1, groups)
        else if (scale.) d1 <- t(scale(t(dset1)))
        if (two.tails) {
            reg <- lapply(1:ncol(d1), function(i, d1, nn, ws) {
                orden <- order(d1[, i])
                tmp <- d1[orden[c(1:nn, (nrow(d1)-nn+1):nrow(d1))], i]
                return(abs(tmp)^ws*sign(tmp))
            }, d1=d1, nn=nn, ws)
        }
        else {
            reg <- lapply(1:ncol(d1), function(i, d1, nn, ws) {
                orden <- order(abs(d1[, i]), decreasing=TRUE)
                tmp <- d1[orden[1:(2*nn)], i]
                return(abs(tmp)^ws)
            }, d1=d1, nn=nn, ws=ws)
        }
        names(reg) <- colnames(d1)
        if (ws==0) {
            d1 <- apply(d1, 2, rank)/(nrow(d1)+1)*2-1
            if (!two.tails) {
                d1 <- abs(d1)*2-1
                d1[d1==(-1)] <- 1-(1/nrow(d1))
            }
            d1 <- qnorm(d1/2+.5)
        }
        genes <- unique(unlist(lapply(reg, names), use.names=FALSE))
        d1 <- t(d1[match(genes, rownames(d1)), ])
        reg <- sapply(reg, function(x, genes) x[match(genes, names(x))], genes=genes)
        reg[is.na(reg)] <- 0
        ww <- colSums(abs(reg))
        reg <- t(t(reg)/ww)
        es <- (d1 %*% reg)*sqrt(ww)
        offdiag <- rowMeans(cbind(es[lower.tri(es)], t(es)[lower.tri(es)]))
        es[lower.tri(es)] <- offdiag
        es <- t(es)
        es[lower.tri(es)] <- offdiag
        class(es) <- "signatureDistance"
        return(es)
    }
    d2 <- dset2
    genes <- rownames(d1)[rownames(d1) %in% rownames(d2)]
    if (ws==0 & length(genes)<(nn*3)) stop("Error, insufficient number of common features to perform analysis", call.=FALSE)
    d1 <- filterRowMatrix(d1, match(genes, rownames(d1)))
    d2 <- filterRowMatrix(d2, match(genes, rownames(d2)))
    if (!is.null(groups)) {
        if (length(groups)==2) {
            d1 <- scaleGroups(dset1, groups[[1]])
            d2 <- scaleGroups(dset2, groups[[2]])
        }
        else {
            d1 <- t(scale(t(dset1)))
            d2 <- t(scale(t(dset2)))
            warning("Not using groups...Expecting a 2 elements list for groups when dset2 is defined", call.=FALSE)
        }
    }
    else if (scale.) {
        d1 <- t(scale(t(dset1)))
        d2 <- t(scale(t(dset2)))
    }
    if (two.tails) {
        reg1 <- lapply(1:ncol(d1), function(i, d1, nn, ws) {
            orden <- order(d1[, i])
            tmp <- d1[orden[c(1:nn, (nrow(d1)-nn+1):nrow(d1))], i]
            return(abs(tmp)^ws*sign(tmp))
        }, d1=d1, nn=nn, ws=ws)
        reg2 <- lapply(1:ncol(d2), function(i, d1, nn, ws) {
            orden <- order(d1[, i])
            tmp <- d1[orden[c(1:nn, (nrow(d1)-nn+1):nrow(d1))], i]
            return(abs(tmp)^ws*sign(tmp))
        }, d1=d2, nn=nn, ws=ws)
    }
    else {
        reg1 <- lapply(1:ncol(d1), function(i, d1, nn, ws) {
            orden <- order(d1[, i], decreasing=TRUE)
            tmp <- d1[orden[1:(2*nn)], i]
            return(abs(tmp)^ws)
        }, d1=d1, nn=nn, ws=ws)
        reg2 <- lapply(1:ncol(d2), function(i, d1, nn, ws) {
            orden <- order(d1[, i], decreasing=TRUE)
            tmp <- d1[orden[1:(2*nn)], i]
            return(abs(tmp)^ws)
        }, d1=d2, nn=nn, ws=ws)
    }
    names(reg1) <- colnames(d1)
    names(reg2) <- colnames(d2)
    if (ws==0) {
        d1 <- apply(d1, 2, rank)/(nrow(d1)+1)*2-1
        d2 <- apply(d2, 2, rank)/(nrow(d2)+1)*2-1
        if (!two.tails) {
            d1 <- abs(d1)*2-1
            d1[d1==(-1)] <- 1-(1/nrow(d1))
            d2 <- abs(d2)*2-1
            d2[d2==(-1)] <- 1-(1/nrow(d2))
        }
        d1 <- qnorm(d1/2+.5)
        d2 <- qnorm(d2/2+.5)
    }
    genes <- unique(unlist(lapply(reg1, names), use.names=FALSE))
    d2 <- t(filterRowMatrix(d2, match(genes, rownames(d2))))
    reg1 <- sapply(reg1, function(x, genes) x[match(genes, names(x))], genes=genes)
    reg1[is.na(reg1)] <- 0
    genes <- unique(unlist(lapply(reg2, names), use.names=FALSE))
    d1 <- t(filterRowMatrix(d1, match(genes, rownames(d1))))
    reg2 <- sapply(reg2, function(x, genes) x[match(genes, names(x))], genes=genes)
    reg2[is.na(reg2)] <- 0
    ww1 <- colSums(abs(reg1))
    ww2 <- colSums(abs(reg2))
    reg1 <- t(t(reg1)/ww1)
    reg2 <- t(t(reg2)/ww2)
    es <- ((d1 %*% reg2)*sqrt(ww1)+t(d2 %*% reg1)*sqrt(ww2))/2
    class(es) <- "signatureDistance"
    return(es)
}

#' Signatures with grouping variable
#' 
#' This function compute signatures using groups information
#' 
#' @param x Numerical matrix with genes in rows and samples in columns
#' @param groups Vector of same length as columns has the dset containing the labels for grouping the samples
#' @return Numeric matrix of signatures (z-scores) with genes in rows and groups in columns
#' @description \code{scaleGroups} compares each group vs. the remaining groups using a Student's t-test
#' @export
#' @examples
#' data(bcellViper, package="bcellViper")
#' res <- scaleGroups(exprs(dset)[, 1:20], rep(1:4, rep(5, 4)))
#' res[1:5, ]
scaleGroups <- function(x, groups) {
    colnames(x) <- groups
    res <- sapply(unique(groups), function(x, dset) {
        x1 <- dset[, colnames(dset)==x]
        x2 <- dset[, colnames(dset)!=x]
        if (is.null(ncol(x1)) | is.null(ncol(x2))) tmp <- rowTtest(x1-x2, alternative="two.sided")
        else tmp <- rowTtest(x1, x2, alternative="two.sided")
        return((qnorm(tmp$p.value/2, lower.tail=FALSE)*sign(tmp$statistic))[, 1])
    }, dset=x)
    mm <- max(abs(res)[is.finite(res)])
    ipos <- which(is.infinite(res))
    res[ipos] <- mm+runif(length(ipos), 0, 10)
    colnames(res) <- unique(groups)
    return(res)
}

#' Distance matrix from signatureDistance objects
#' 
#' This function transforms a signatureDistance object into a dist object
#' 
#' @param m signatureDistance object
#' @param diag parameter included for compatibility
#' @param upper parameter included for compatibility
#' @return Object of class dist
#' @method as.dist signatureDistance
#' @export
#' @importFrom stats as.dist
as.dist.signatureDistance <- function(m, diag=FALSE, upper=FALSE) {
    m <- scale(m)
    class(m) <- "matrix"
    return(as.dist(1-m))
}

#' Scaling of signatureDistance objects
#' 
#' This function scales the signatureDistance so its range is (-1, 1)
#' 
#' @param x signatureDistance object
#' @param center Not used, given for compatibility with the generic funtion scale
#' @param scale Not used, given for compatibility with the generic function scale
#' @return Scaled signatureDistance object
#' @method scale signatureDistance
#' @export
scale.signatureDistance <- function(x, center=TRUE, scale=TRUE) {
    if (ncol(x) != nrow(x)) return(x/max(abs(x)))
    mm <- (matrix(diag(x), nrow(x), ncol(x)) + matrix(diag(x), nrow(x), ncol(x), byrow=TRUE))/2
    return(x/mm)
}

#' VIPER similarity
#' 
#' This function computes the similarity between VIPER signatures
#' 
#' @param x Numeric matrix containing the VIPER results with samples in columns and regulators in rows
#' @param nn Optional number of top regulators to consider for computing the similarity
#' @param ws Number indicating the weighting exponent for the signature, or vector of 2 numbers indicating the inflection point and the value corresponding to a weighting score of .1 for a sigmoid transformation, only used if \code{nn} is ommited
#' @param method Character string indicating whether the most active (greater), less active (less) or both tails (two.sided) of the signature should be used for computing the similarity
#' @description If ws is a single number, weighting is performed using an exponential function. If ws is a 2 numbers vector, weighting is performed with a symmetric sigmoid function using the first element as inflection point and the second as trend.
#' @return signatureDistance object
#' @examples
#' data(bcellViper, package="bcellViper")
#' dd <- viperSimilarity(exprs(dset))
#' dd[1:5, 1:5]
#' scale(dd)[1:5, 1:5]
#' as.matrix(as.dist(dd))[1:5, 1:5]
#' @export

viperSimilarity <- function(x, nn=NULL, ws=c(4, 2), method=c("two.sided", "greater", "less")) {
    method <- match.arg(method)
    if (min(x, na.rm=TRUE)>=0) {
        x <- apply(x, 2, function(x) rank(x, na.last="keep"))
        x <- qnorm(t(t(x)/(colSums(!is.na(x))+1)))
    }
    x[is.na(x)] <- 0
    xw <- x
    if (is.null(nn)) {
        if (length(ws)==1) {
            switch(method,
                   greater={
                        xw[xw<0] <- 0
                        xw <- t(t(xw)/apply(x, 2, max))
                   },
                   less={
                        xw[xw>0] <- 0
                        xw <- t(t(xw)/apply(abs(x), 2, max))
                   },
                   two.sided={
                        xw <- t(t(xw)/apply(abs(x), 2, max))
                   })
            xw <- sign(xw)*abs(xw)^ws
        }
        else {
            ws[2] <- 1/(ws[2]-ws[1])*log(1/.9-1)
            switch(method,
                   greater={
                       xw[xw<0] <- 0
                   },
                   less={
                       xw[xw>0] <- 0
                   },
                   two.sided={
                   })
            xw <- sign(xw)*sigT(abs(xw), ws[2], ws[1])
        }
    }
    else {
        switch(method,
               greater={
                   xw <- apply(xw, 2, function(x, nn) {
                        pos <- rank(-x, ties.method="random")
                        x[pos>nn] <- NA
                        return(x)
                   }, nn=nn)
               },
               less={
                   xw <- apply(xw, 2, function(x, nn) {
                        pos <- rank(x, ties.method="random")
                        x[pos>nn] <- NA
                        return(x)
                   }, nn=nn)
               },
               two.sided={
                   xw <- apply(xw, 2, function(x, nn) {
                        pos <- rank(x, ties.method="random")
                        x[pos>nn & pos < (length(x)-nn+1)] <- NA
                        return(x)
                   }, nn=round(nn/2))
               })
        xw <- x/abs(xw)
        xw[is.na(xw)] <- 0
    }
    nes <- sqrt(colSums(xw^2))
    xw <- scale(xw, center=FALSE, scale=colSums(abs(xw)))
    t2 <- qnorm(apply(x, 2, rank)/(nrow(x)+1))
    vp <- t(xw) %*% t2
    vp <- vp * nes
    tmp <- cbind(vp[lower.tri(vp)], t(vp)[lower.tri(vp)])
    tmp <- rowSums(tmp*tmp^2)/rowSums(tmp^2)
    vp[lower.tri(vp)] <- tmp
    vp <- t(vp)
    vp[lower.tri(vp)] <- tmp
    class(vp) <- "signatureDistance"
    return(vp)
}
