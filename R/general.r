##########
#' Coefficient of variation filter
#'
#' This function filter redundant probes based on the highest coefficient of variation
#'
#' @param expset Expression set or Matrix containing the gene expression data, with samples in columns and probes in rows. The \code{colnames} attribute should contain the sample names and the \code{rownames} attribute should contain the unique geneIDs
#' @param ... Additional parameters added to keep compatibility
#' @return CV filtered dataset
#' @examples
#' data(bcellViper, package="bcellViper")
#' d1 <- exprs(dset)
#' tmp <- rownames(d1)
#' tmp[round(runif(10, 1, length(tmp)))] <- tmp[1]
#' rownames(d1) <- tmp
#' dim(d1)
#' d1 <- filterCV(d1)
#' dim(d1)
#' @export
#' @docType methods
#' @rdname filterCV-methods
setGeneric("filterCV", function(expset, ...) standardGeneric("filterCV"))

#' @rdname filterCV-methods
#' @aliases filterCV,matrix-method
setMethod("filterCV", "matrix", function (expset) {
    repet <- tapply(rownames(expset), factor(rownames(expset)), length)
    if (max(repet)>1) {
        d <- expset[rownames(expset) %in% (names(repet)[repet==1]),]  
        d1 <- expset[rownames(expset) %in% (names(repet)[repet>1]),]
        d2 <- tapply(1:nrow(d1),factor(rownames(d1)), function(pos, d1, cv) {
            list(name=rownames(d1)[pos][which.max(cv[pos])],value=d1[pos,][which.max(cv[pos]),])
        }, d1=d1, cv=frcv(d1))
        d3 <- t(sapply(d2, function(x) x$value))
        if (nrow(d3) != length(d2)) d3 <- t(d3)
        rownames(d3) <- sapply(d2, function(x) x$name)
        expset <- rbind(d,d3)
    }
    expset
})

#' @rdname filterCV-methods
#' @aliases filterCV,ExpressionSet-method
setMethod("filterCV", "ExpressionSet", function(expset) {
    exprs(expset) <- filterCV(exprs(expset))
    return(expset)
})

##########
#' Null model by sample permutation testing
#'
#' This function performs sample permutation and t-test to generate a null model
#'
#' @param x ExpressionSet object or Matrix containing the test dataset
#' @param ... Additional parameters added to keep compatibility
#' @return Matrix of z-scores with genes in rows and permutations in columns
#' @seealso \code{\link{msviper}}, \code{\link{viper}}
#' @export
#' @docType methods
#' @rdname ttestNull-methods
setGeneric("ttestNull", function(x, ...) standardGeneric("ttestNull"))

#' @param y Matrix containing the reference dataset
#' @param per Integer indicating the number of permutations
#' @param repos Logical, whether the permutations should be performed with reposition
#' @param seed Integer indicating the seed for the permutations, 0 for disable it
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @examples
#' data(bcellViper, package="bcellViper")
#' d1 <- exprs(dset)
#' dnull <- ttestNull(d1[, 1:10], d1[, 11:20], per=100)
#' dim(dnull)
#' plot(density(dnull))
#' @rdname ttestNull-methods
#' @aliases ttestNull,matrix-method
setMethod("ttestNull", c(x="matrix"), function(x, y, per=1000, repos=TRUE, seed=1, cores=1, verbose=TRUE) {
    if (seed>0) set.seed(round(seed))
    pb <- NULL
    if (verbose) {
        message(date(), "\nComputing the null model distribution by ", per, " permutations.")
    }
    if (cores>1) {
        res <- .viperParallelApply(1:per, function(i, x, y, repos) {
            expset <- cbind(x, y)
            repeat{
                sorder <- sample(ncol(expset), replace=repos)
                if (length(unique(sorder[1:ncol(x)]))>1 & length(unique(sorder[-(1:ncol(x))]))>1) break
                if (verbose) message("-", appendLF=FALSE)
            }
            x1 <- filterColMatrix(expset, sorder[1:ncol(x)])
            y1 <- filterColMatrix(expset, sorder[-(1:ncol(x))])
            largo <- rowSums(!is.na(x1))
            largoy <- rowSums(!is.na(y1))
            t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
            t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
            names(t) <- rownames(x)
            return(t)
        }, x=x, y=y, repos=repos, mc.cores=cores)
        res <- sapply(res, function(x) x)
    }
    else {
        if (verbose) pb <- txtProgressBar(max=per, style=3)
        res <- sapply(1:per, function(i, x, y, repos, pb, verbose) {
            if (verbose) setTxtProgressBar(pb, i)
            expset <- cbind(x, y)
            repeat{
                sorder <- sample(ncol(expset), replace=repos)
                if (length(unique(sorder[1:ncol(x)]))>1 & length(unique(sorder[-(1:ncol(x))]))>1) break
                if (verbose) message("-", appendLF=FALSE)
            }
            x1 <- filterColMatrix(expset, sorder[1:ncol(x)])
            y1 <- filterColMatrix(expset, sorder[-(1:ncol(x))])
            largo <- rowSums(!is.na(x1))
            largoy <- rowSums(!is.na(y1))
            t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
            t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
            names(t) <- rownames(x)
            return(t)
        }, x=x, y=y, repos=repos, pb=pb, verbose=verbose)
    }
    colnames(res) <- 1:per
    if (verbose) message("\n", date())
    return(res)
})

#' @param pheno Character string indicating the phenotype data to use
#' @param group1 Vector of character strings indicating the category from phenotype \code{pheno} to use as test group
#' @param group2 Vector of character strings indicating the category from phenotype \code{pheno} to use as control group
#' @examples
#' data(bcellViper, package="bcellViper")
#' dnull <- ttestNull(dset, "description", "CB", "CC", per=100)
#' dim(dnull)
#' plot(density(dnull))
#' @rdname ttestNull-methods
#' @aliases ttestNull,ExpressionSet-method
setMethod("ttestNull", c(x="ExpressionSet"), function(x, pheno, group1, group2, per=1000, repos=TRUE, seed=1, verbose=TRUE) {
    pos1 <- pData(x)[[pheno]] %in% group1
    pos2 <- pData(x)[[pheno]] %in% group2
    if (length(pos1)==0) stop(paste(pheno, " is not present in the ExpressionSet object", sep=""), call.=FALSE)
    if (length(which(pos1))==0) stop(paste(group1, " was nor found in ", pheno, sep=""), call.=FALSE)
    if (length(which(pos2))==0) stop(paste(group2, " was nor found in ", pheno, sep=""), call.=FALSE)
    ttestNull(exprs(x)[, pos1], exprs(x)[, pos2], per=per, repos=repos, seed=seed, verbose=verbose)
})

##########
#' Bootstrapped signature by t-test
#'
#' This function generates a bootstrapped signature matrix by t-test
#'
#' @param x Matrix containing the test dataset
#' @param y Matrix containing the reference dataset
#' @param per Integer indicating the number of permutations
#' @param seed Integer indicating the seed for the permutations, 0 for disable it
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical whether progress should be reported
#' @param ... Additional parameters added to keep compatibility
#' @return Matrix of z-scores with genes in rows and permutations in columns
#' @seealso \code{\link{msviper}}
#' @export
#' @docType methods
#' @rdname bootstrapTtest-methods
setGeneric("bootstrapTtest", function(x, ...) standardGeneric("bootstrapTtest"))

#' @examples
#' data(bcellViper, package="bcellViper")
#' d1 <- exprs(dset)
#' sig <- bootstrapTtest(d1[, 1:10], d1[, 11:20], per=100)
#' dim(sig)
#' plot(density(sig[1907, ]))
#' @rdname bootstrapTtest-methods
#' @aliases bootstrapTtest,matrix-method
setMethod("bootstrapTtest", c(x="matrix"), function(x, y, per=100, seed=1, cores=1, verbose=TRUE) {
    if (seed>0) set.seed(round(seed))
    pb <- NULL
    if (verbose) {
        message(date(), "\nComputing the bootstrapped signatures by ", per, " permutations.")
    }
    if (cores>1) {
        res <- .viperParallelApply(1:per, function(i, x, y) {
            x1 <- filterColMatrix(x, sample(ncol(x), replace=TRUE))
            y1 <- filterColMatrix(y, sample(ncol(y), replace=TRUE))
            largo <- rowSums(!is.na(x1))
            largoy <- rowSums(!is.na(y1))
            t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
            t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
            names(t) <- rownames(x1)
            return(t)
        }, x=x, y=y, mc.cores=cores)
        res <- sapply(res, function(x) x)
    }
    else {
        if (verbose) pb <- txtProgressBar(max=per, style=3)
        res <- sapply(1:per, function(i, x, y, pb, verbose) {
            if (verbose) setTxtProgressBar(pb, i)
            x1 <- filterColMatrix(x, sample(ncol(x), replace=TRUE))
            y1 <- filterColMatrix(y, sample(ncol(y), replace=TRUE))
            largo <- rowSums(!is.na(x1))
            largoy <- rowSums(!is.na(y1))
            t <- ((rowMeans(x1, na.rm=TRUE) - rowMeans(y1, na.rm=TRUE))/sqrt(((largo - 1) *  rowVars(x1) + (largoy - 1) * rowVars(y1))/(largo + largoy - 2))/sqrt(1/largo + 1/largoy))[, 1]
            t <- qnorm(pt(abs(t), largo + largoy - 2, lower.tail = FALSE), lower.tail=FALSE)*sign(t)
            names(t) <- rownames(x1)
            return(t)
        }, x=x, y=y, pb=pb, verbose=verbose)
    }
    colnames(res) <- 1:per
    if (verbose) message("\n", date())
    return(res)
})

#' @param pheno Character string indicating the phenotype data to use
#' @param group1 Vector of character strings indicating the category from phenotype \code{pheno} to use as test group
#' @param group2 Vector of character strings indicating the category from phenotype \code{pheno} to use as control group
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- bootstrapTtest(dset, "description", "CB", "N", per=100)
#' dim(sig)
#' plot(density(sig[1907, ]))
#' @rdname bootstrapTtest-methods
#' @aliases bootstrapTtest,ExpressionSet-method
setMethod("bootstrapTtest", c(x="ExpressionSet"), function(x, pheno, group1, group2, per=100, seed=1, verbose=TRUE) {
    pos1 <- pData(x)[[pheno]] %in% group1
    pos2 <- pData(x)[[pheno]] %in% group2
    if (length(pos1)==0) stop(paste(pheno, " is not present in the ExpressionSet object", sep=""), call.=FALSE)
    if (length(which(pos1))==0) stop(paste(group1, " was nor found in ", pheno, sep=""), call.=FALSE)
    if (length(which(pos2))==0) stop(paste(group2, " was nor found in ", pheno, sep=""), call.=FALSE)
    bootstrapTtest(exprs(x)[, pos1], exprs(x)[, pos2], per=per, seed=seed, verbose=verbose)
})

##########
#' msVIPER annotation change
#'
#' This function changes the annotation of genes in msviper objects
#'
#' @param mobj msviper object generated by \code{msviper} function
#' @param annot Vector os character strings containing the gene names and gene identifiers as vector names attribute
#' @param complete Logical, whether the signature and target names should be also transformed
#' @return msviper object with updated annotations
#' @seealso \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- rowTtest(dset, "description", "CB", "N")$statistic
#' mra <- msviper(sig, regulon)
#' tmp <- unique(c(names(mra$regulon), rownames(mra$signature)))
#' annot <- 1:length(tmp)
#' names(annot) <- tmp
#' plot(mra, cex=.7)
#' mra <- msviperAnnot(mra, annot)
#' plot(mra, cex=.7)
#' @export

msviperAnnot <- function(mobj, annot, complete=TRUE) {
    names(mobj$regulon) <- comNames(names(mobj$regulon), annot)
    mobj$regulon <- mobj$regulon[!is.na(names(mobj$regulon))]
    mobj$es <- lapply(mobj$es[sapply(mobj$es, length)>0], function(x, annot) {
        if (is.null(ncol(x))) {
            names(x) <- comNames(names(x), annot)
            x <- x[!is.na(names(x))]
        }
        else {
            rownames(x) <- comNames(rownames(x), annot)
            x <- x[!is.na(rownames(x)), ]
        }
        return(x)
    }, annot=annot)
    if (complete) {
        if (is.null(dim(mobj$signature))) {
            names(mobj$signature) <- annot[match(names(mobj$signature), names(annot))]
            mobj$signature <- mobj$signature[!is.na(names(mobj$signature))]
        }
        else {
            rownames(mobj$signature) <- annot[match(rownames(mobj$signature), names(annot))]
            mobj$signature <- filterRowMatrix(mobj$signature, !is.na(rownames(mobj$signature)))
        }
        mobj$regulon <- lapply(mobj$regulon, function(x, annot) {
            names(x$tfmode) <- annot[match(names(x$tfmode), names(annot))]
            filtro <- !is.na(names(x$tfmode))
            x$tfmode <- x$tfmode[filtro]
            x$likelihood <- x$likelihood[filtro]
            return(x)
        }, annot=annot)
        if (!is.null(mobj$nullmodel)) {
            rownames(mobj$nullmodel) <- annot[match(rownames(mobj$nullmodel), names(annot))]
            mobj$nullmodel <- mobj$nullmodel[!is.na(rownames(mobj$nullmodel)), ]
        }
    }
    return(mobj)
}

##########
#' Combinatorial annotation
#'
#' This function convers combinatorial annotations
#'
#' @param x Character vector of gene name combinations, where the combinations are separated by --
#' @param annot Vector of gene names with geneID as \code{names} attribute
#' @return Converted annotations
#' @seealso \code{\link{msviper}}

comNames <- function(x, annot) {
    sapply(strsplit(as.character(x), "--"), function(x, annot) {
        paste(annot[match(x, names(annot))], collapse="--")
    }, annot=annot)
}

rowVars <- function(x) {
    ave <- rowMeans(x, na.rm=TRUE)
    pos <- which(is.na(x))
    largo <- rowSums(!is.na(x))
    x[pos] <- rep(ave, ncol(x))[pos]
    (x - ave)^2 %*% rep(1, ncol(x))/(largo - 1)
}

##########
#' Prune Regulons
#' 
#' This function limits the maximum size of the regulons
#' 
#' @param regulon Object of class regulon
#' @param cutoff Number indicating the maximum size for the regulons (maximum number of target genes)
#' @param adaptive Logical, whether adaptive size should be used (i.e. sum(likelihood^2))
#' @param eliminate Logical whether regulons smalles than \code{cutoff} should be eliminated
#' @param wm Optional numeric vector of weights (0; 1) for the genes
#' @return Prunned regulon
#' @seealso \code{\link{viper}}, \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' hist(sapply(regulon, function(x) sum(x$likelihood)/max(x$likelihood)), nclass=20)
#' preg <- pruneRegulon(regulon, 400)
#' hist(sapply(preg, function(x) sum(x$likelihood)/max(x$likelihood)), nclass=20)
#' @export
pruneRegulon <- function(regulon, cutoff=50, adaptive=TRUE, eliminate=FALSE, wm=NULL) {
    if (adaptive) {
        regulon <- lapply(regulon, function(x, cutoff, wm) {
            likelihood <- x$likelihood
            if (!is.null(wm)) {
                wm <- wm[match(names(x$tfmode), names(wm))]
                wm[is.na(wm)] <- 0
                likelihood <- likelihood * wm
            }
            pos <- order(likelihood, decreasing=TRUE)
            ws <- (likelihood/max(likelihood))^2
            pos <- pos[cumsum(ws[pos])<=cutoff]
            return(list(tfmode=x$tfmode[pos], likelihood=x$likelihood[pos]))
        }, cutoff=cutoff, wm=wm)
    }
    else {
        regulon <- lapply(regulon, function(x, cutoff) {
            pos <- order(x$likelihood, decreasing=TRUE)
            pos <- pos[1:min(length(pos), cutoff)]
            return(list(tfmode=x$tfmode[pos], likelihood=x$likelihood[pos]))
        }, cutoff=cutoff)
        if (eliminate) regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=cutoff]
    }
    class(regulon) <- "regulon"
    return(regulon)
}

##########
#' Integrate signatures
#' 
#' This function integrates signatures represented as columns in the input matrix using self-weighting average
#' 
#' @param signature Numeric matrix containing the signatures as z-scores or NES, genes in rows and signatures in columns
#' @param score Number indicating the exponent score for the weight
#' @return Vector containing the integrated signatures
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- bootstrapTtest(dset, "description", "CB", "N", per=100)
#' isig <- integrateSignatures(sig)
#' plot(density(sig))
#' lines(density(isig, adj=1.5), col="red")
#' @export
integrateSignatures <- function(signature, score=1) {
    w <- abs(signature)^score
    w <- w/rowSums(w)
    return(rowSums(signature*w))
}

##########
#' Variance of rows for arrays with NA values
#'
#' This function computes the variance by rows ignoring NA values
#'
#' @param x Numeric matrix
#' @return 1-column matrix with the variance by row results
#' @examples
#' data(bcellViper, package="bcellViper")
#' tmp <- exprs(dset)[1:10, ]
#' tmp[round(runif(100, 1, length(tmp)))] <- NA
#' frvarna(tmp)
#' @export
frvarna <- function(x) {
    # Missing values imputed with the row mean contribute zero to the
    # numerator, so this is equivalent to sample variance over observed
    # values. matrixStats computes each row in compiled code without copying
    # and imputing the full input matrix.
    out <- matrix(matrixStats::rowVars(x, na.rm=TRUE), ncol=1L,
                  dimnames=list(rownames(x), NULL))
    out
}

# Compiled row ranks with base::rank()'s default treatment of missing values:
# missing entries are tied at the end rather than left as NA.
fastRowRanks <- function(x) {
    out <- matrixStats::rowRanks(x, ties.method="average")
    missing <- is.na(x)
    if (any(missing)) {
        idx <- which(missing, arr.ind=TRUE)
        nmissing <- rowSums(missing)
        missingOrder <- matrixStats::rowCumsums(missing)
        out[idx] <- ncol(x) - nmissing[idx[, 1L]] +
            missingOrder[idx]
    }
    out
}

##########
#' Variance of columns for arrays with NA values
#'
#' This function computes the variance by columns ignoring NA values
#'
#' @param x Numeric matrix
#' @return 1-column matrix with the variance by column results
#' @examples
#' data(bcellViper, package="bcellViper")
#' tmp <- exprs(dset)[, 1:10]
#' tmp[round(runif(100, 1, length(tmp)))] <- NA
#' fcvarna(tmp)
#' @export
fcvarna <- function(x) frvarna(t(x))

##########
#' Student's t-test for rows
#' 
#' This function performs a Student's t-test on each row of a matrix
#' 
#' @param x ExpressionSet object or Numerical matrix containing the test samples
#' @param ... Additional parameters added to keep compatibility
#' @return List of Student-t-statistic (\code{statistic}) and p-values (\code{p.value})
#' @export
#' @docType methods
#' @rdname rowTtest-methods
setGeneric("rowTtest", function(x, ...) standardGeneric("rowTtest"))

#' @param y Optional numerical matrix containing the reference samples. If ommited \code{x} will be tested against mean = \code{mu}
#' @param mu Number indicating the alternative hypothesis when \code{y} is ommited
#' @param alternative Character string indicating the tail for the test, either two.sided, greater or lower
#' @examples
#' data(bcellViper, package="bcellViper")
#' d1 <- exprs(dset)
#' res <- rowTtest(d1[, 1:10], d1[, 11:20])
#' res$statistic[1:5, ]
#' res$p.value[1:5, ]
#' @rdname rowTtest-methods
#' @aliases rowTtest,matrix-method
setMethod("rowTtest", c(x="matrix"), function(x, y=NULL, mu=0, alternative="two.sided") {
    largo <- rowSums(!is.na(x))
    if (is.null(y)) {
        t <- (rowMeans(x, na.rm=TRUE)-mu)/sqrt(frvarna(x)/largo)
        pval <- switch(pmatch(alternative, c("two.sided", "greater", "less")),
        pt(abs(t),largo-1,lower.tail=FALSE)*2,
        pt(t, largo-1, lower.tail=FALSE),
        pt(t, largo-1, lower.tail=TRUE))
        list(statistic=t, p.value=pval)
    }
    else {
        largoy <- rowSums(!is.na(y))
        t <- (rowMeans(x, na.rm=TRUE)-rowMeans(y, na.rm=TRUE))/sqrt(((largo-1)*frvarna(x)+(largoy-1)*frvarna(y))/(largo+largoy-2))/sqrt(1/largo+1/largoy)
        pval <- switch(pmatch(alternative, c("two.sided", "greater", "less")),
        pt(abs(t),largo+largoy-2,lower.tail=FALSE)*2,
        pt(t, largo+largoy-2, lower.tail=FALSE),
        pt(t, largo+largoy-2, lower.tail=TRUE))
        list(statistic=t, p.value=pval)
    }
})

#' @param pheno Character string indicating the phenotype data to use
#' @param group1 Vector of character strings indicating the category from phenotype \code{pheno} to use as test group
#' @param group2 Vector of character strings indicating the category from phenotype \code{pheno} to use as control group
#' @examples
#' data(bcellViper, package="bcellViper")
#' res <- rowTtest(dset, "description", "CB", "N")
#' res$statistic[1:5, ]
#' res$p.value[1:5, ]
#' @rdname rowTtest-methods
#' @aliases rowTtest,ExpressionSet-method
setMethod("rowTtest", c(x="ExpressionSet"), function(x, pheno, group1, group2=NULL, mu=0, alternative="two.sided") {
    if (is.null(group2)) {
        pos <- pData(x)[[pheno]] %in% group1
        if (length(pos)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
        if (length(which(pos))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
        return(rowTtest(exprs(x)[, pos], mu=mu, alternative=alternative))
    }
    pos1 <- pData(x)[[pheno]] %in% group1
    pos2 <- pData(x)[[pheno]] %in% group2
    if (length(pos1)==0) stop(paste(pheno, " was not found in the ExpressionSet Object", sep=""), call.=FALSE)
    if (length(which(pos1))==0) stop(paste(group1, " was not found in ", pheno, sep=""), call.=FALSE)
    if (length(which(pos2))==0) stop(paste(group2, " was not found in ", pheno, sep=""), call.=FALSE)
    return(rowTtest(exprs(x)[, pos1], exprs(x)[, pos2], alternative=alternative))
})

##########
#' Coeficient of variations for rows
#'
#' This function computes the coefficient of variation (CV) by rows
#'
#' @param x Numeric matrix
#' @return 1-column matrix with the coefficient of variation by row results
#' @examples
#' data(bcellViper, package="bcellViper")
#' tmp <- exprs(dset)[1:10, ]
#' tmp[round(runif(100, 1, length(tmp)))] <- NA
#' frcv(tmp)
#' @export

frcv <- function(x) sqrt(frvarna(x))/rowMeans(x, na.rm=TRUE)

##########
#' Loading expression sets
#' 
#' This function load an expression file into a matrix
#' 
#' @param filename Character string indicating the name of the expression file
#' @return List containing a numeric matrix of expression data with samples in columns and probes in rows; and a vector of gene mapping annotations
loadExpset <- function(filename) {
    tmp <- strsplit(readLines(filename), "\t")
    d1 <- t(sapply(tmp[-1], function(x) as.numeric(x[-(1:2)])))
    colnames(d1) <- tmp[[1]][-(1:2)]
    rownames(d1) <- sapply(tmp[-1], function(x) x[1])
    annot <- sapply(tmp[-1], function(x) x[2])
    names(annot) <- rownames(d1)
    return(list(expset=d1, annot=annot))
}

#' msVIPER class
#' 
#' This function generates an instance of the msviper class from a signature, NES signature and regulon object
#' 
#' @param nes Numeric vector of NES values
#' @param signature Numeric vector of gene expression signature
#' @param regulon Instance of class regulon
#' @param nullmodel Optional matrix containing the signatures for the null model 
#' @return msviper class object
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- rowTtest(dset, "description", c("CB", "CC"), "N")$statistic
#' mra <- msviper(sig, regulon)
#' mra1 <- msviperClass(mra$es$nes, sig, regulon)
#' summary(mra1)
#' plot(mra1)
#' @export
msviperClass <- function(nes, signature, regulon, nullmodel=NULL) {
    if (!is.null(ncol(nes))) {
        warning("nes is a matrix, only the first column will be used", call.=FALSE)
        nes <- nes[, 1]
    }
    if (!is.null(ncol(signature))) {
        warning("signature is a matrix, only the first column will be used", call.=FALSE)
        signature <- signature[, 1]
    }
    genes <- intersect(names(nes), names(regulon))
    regulon <- regulon[match(genes, names(regulon))]
    regulon <- lapply(regulon, function(x, genes) {
        pos <- which(names(x$tfmode) %in% genes)
        list(tfmode=x$tfmode[pos], likelihood=x$likelihood[pos])
    }, genes=names(signature))
    class(regulon) <- "regulon"
    nes <- nes[match(genes, names(nes))]
    res <- list(signature=matrix(signature, length(signature), 1, dimnames=list(names(signature), 1)), regulon=regulon, es=list(nes=nes, nes.se=NULL, size=sapply(regulon, function(x) length(x$tfmode)), p.value=pnorm(abs(nes), lower.tail=FALSE)*2, nes.bt=matrix(nes, length(nes), 1, dimnames=list(names(nes), 1))), param=list(minsize=25, adaptive.size=FALSE), nullmodel=nullmodel)
    class(res) <- "msviper"
    return(res)
}
