# Functions for computing RPT-activity

#' viperRPT
#' 
#' This function computes residual post-translational activity
#' 
#' @param vipermat Numeric matrix containing the viper protein activity inferences
#' @param expmat Numeric matrix or expressionSet containing the expression data
#' @param weights List of numeric matrix of sample weights
#' @param method Character string indicating the method to use, either rank, lineal or spline
#' @param robust Logical, whether the contribution of outliers is down-weighted by using a gaussian kernel estimate for the join probability density
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @return Matrix of RPT-activity values
#' @seealso \code{\link{viper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' vipermat <- viper(dset, regulon)
#' rpt <- viperRPT(vipermat, dset)
#' rpt[1:5, 1:5]
#' @export
viperRPT <- function(vipermat, expmat, weights=matrix(1, nrow(vipermat), ncol(vipermat), dimnames=list(rownames(vipermat), colnames(vipermat))), method=c("spline", "lineal", "rank"), robust=FALSE, cores=1) {
    method <- match.arg(method)
    if (is(vipermat, "ExpressionSet")) vipermat <- exprs(vipermat)
    if (is(expmat, "ExpressionSet")) expmat <- exprs(expmat)
    if (is(weights, "list")) {
        weights <- t(sapply(weights, function(x, samp) {
            as.numeric(samp %in% x)
        }, samp=colnames(vipermat)))
    }
# Compaibilizing the matrixes
    genes <- intersect(rownames(vipermat), intersect(rownames(expmat), rownames(weights)))
    samp <- intersect(colnames(vipermat), intersect(colnames(expmat), colnames(weights)))
    vipermat <- filterColMatrix(filterRowMatrix(vipermat, match(genes, rownames(vipermat))), match(samp, colnames(vipermat)))
    expmat <- filterColMatrix(filterRowMatrix(expmat, match(genes, rownames(expmat))), match(samp, colnames(expmat)))
    weights <- filterColMatrix(filterRowMatrix(weights, match(genes, rownames(weights))), match(samp, colnames(weights)))
    # For ordinary complete data, avoid fitting one lm() object per gene.
    # Weighted simple regression with an intercept has a closed form that can
    # be evaluated over all rows with compiled matrix row operations.
    if (!robust && method %in% c("rank", "lineal") &&
        all(is.finite(vipermat)) && all(is.finite(expmat)) &&
        all(is.finite(weights)) && all(weights >= 0)) {
        x <- if (method == "rank") fastRowRanks(expmat) else expmat
        y <- if (method == "rank") fastRowRanks(vipermat) else vipermat
        sw <- rowSums(weights)
        xbar <- rowSums(weights * x) / sw
        ybar <- rowSums(weights * y) / sw
        dx <- x - xbar
        dy <- y - ybar
        sxx <- rowSums(weights * dx^2)
        if (all(sw > 0 & sxx > 0)) {
            slope <- rowSums(weights * dx * dy) / sxx
            intercept <- ybar - slope * xbar
            res <- y - intercept - slope * x
            dimnames(res) <- dimnames(vipermat)
            return(res)
        }
    }
# fitting the models
    switch(method,
    rank={
        res <- .viperParallelApply(1:nrow(vipermat), function(i, vipermat, expmat, weights, robust) {
            x <- rank(expmat[i, ])
            y <- rank(vipermat[i, ])
            w <- weights[i, ]
            if (robust) {
                tmp <- approxk2d(cbind(x, y))
                w <- w*(tmp/max(tmp))^2
            }
            residuals(lm(y~x, data=list(x=x, y=y), weights=w))
        }, vipermat=vipermat, expmat=expmat, weights=weights, robust=robust, mc.cores=cores)
    }, 
    lineal={
        res <- .viperParallelApply(1:nrow(vipermat), function(i, vipermat, expmat, weights, robust) {
            x <- expmat[i, ]
            y <- vipermat[i, ]
            w <- weights[i, ]
            if (robust) {
                tmp <- approxk2d(cbind(x, y))
                w <- w*(tmp/max(tmp))^2
            }
            residuals(lm(y~x, data=list(x=x, y=y), weights=w))
        }, vipermat=vipermat, expmat=expmat, weights=weights, robust=robust, mc.cores=cores)    
    },
    spline={
        res <- .viperParallelApply(1:nrow(vipermat), function(i, vipermat, expmat, weights, robust) {
            x <- expmat[i, ]
            y <- vipermat[i, ]
            w <- weights[i, ]
            if (robust) {
                tmp <- approxk2d(cbind(x, y))
                w <- w*(tmp/max(tmp))^2
            }
            residuals(smooth.spline(x=x, y=y, w=w, spar=1.2, all.knots=TRUE))
        }, vipermat=vipermat, expmat=expmat, weights=weights, robust=robust, mc.cores=cores)        
    })
    names(res) <- rownames(vipermat)
    res <- t(sapply(res, function(x) x))
    colnames(res) <- colnames(vipermat)
    return(res)
}

#' approxk2d
#' 
#' This function uses a gaussian kernel to estimate the joint density distribution at the specified points
#' 
#' @param x Matrix of x and y points
#' @param gridsize number or vector indicating the size of the greed where to estimate the density
#' @param pos Matrix of coordinates to evaluate the density
#' @return Vector of density estimates
#' @examples
#' x <- rnorm(500)
#' y <- x+rnorm(500)
#' kde2 <- approxk2d(cbind(x, y))
#' plot(x, y, pch=20, col=hsv(0, kde2/max(kde2), 1))
#' @export

approxk2d <- function(x, gridsize=128, pos=x) {
    bw <- diff(apply(x, 2, quantile, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))/25
    bw[bw == 0] <- 1
    if (length(gridsize)==1) gridsize <- rep(gridsize, 2)
    tmp <- bkde2D(x, bandwidth=bw, gridsize=gridsize[1:2])
    return(e1071::interpolate(pos, tmp$fhat, adims=list(tmp$x1, tmp$x2)))
}
    
