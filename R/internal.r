#Internal functions for viper package
#Check and update regulons to version 2
updateRegulon <- function(regul) {
    if (is.null(names(regul[[1]]))) {
        tmp <- lapply(regul, function(x) {
            tmp <- rep(0, length(x))
            names(tmp) <- x
            list(tfmode=tmp, likelihood=rep(1, length(tmp)))
        })
        return(tmp)
    }
    if (names(regul[[1]])[1]=="tfmode") return(regul)
    return(lapply(regul, function(x) list(tfmode=x, likelihood=rep(1, length(x)))))
}

#' Proportionally Weighted Enrichment Analysis for gene-set groups
#'
#' This function performs a Proportionally Weighted Enrichment Analysis on groups of gene-sets
#'
#' @param rlist Named vector containing the scores to rank the expression profile or matrix where columns contains bootstraped signatures
#' @param groups List of gene-sets (regulons), each component is a list of two vectors: \emph{TFmode} containing the TFMoA index (-1; 1) and \emph{likelihood} containing the interaction relative likelihood
#' @param nullpw Numerical matrix representing the null model, with genes as rows (geneID as rownames) and permutations as columns
#' @param alternative Character string indicating the alternative hypothesis, either two.sided, greater or less
#' @param per Integer indicating the number of permutations for the genes in case "nullpw" is ommited
#' @param minsize Integer indicating the minimum size for the regulons
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A list containing four matrices:
#' \describe{
#' \item{es}{Enrichment score}
#' \item{nes}{Normalized Enrichment Score}
#' \item{size}{Regulon size}
#' \item{p.value}{Enrichment p.value}
#' }

groupPwea3 <- function(rlist, groups, nullpw=NULL, alternative=c("two.sided", "less", "greater"), per=0, minsize=5, cores=1, verbose=TRUE) {
    alternative <- match.arg(alternative)
    if (is.vector(rlist)) rlist <- matrix(rlist, length(rlist), 1, dimnames=list(names(rlist), NULL))
    if (!is.list(groups) | names(groups)[1]=="tfmode") groups <- list(set=groups)
    groups <- updateRegulon(groups)
    rlist2 <- apply(rlist, 2, rank)/(nrow(rlist)+1)*2-1
    rlist1 <- abs(rlist2)*2-1
    rlist1 <- rlist1+(1-max(rlist1))/2
    rlist1 <- qnorm(rlist1/2+.5)
    rlist2 <- qnorm(rlist2/2+.5)
    groups <- lapply(groups, function(x, genes) {
        tmp <- names(x$tfmode) %in% genes
        x$tfmode <- x$tfmode[tmp]
        if(length(x$likelihood)==length(tmp)) x$likelihood <- x$likelihood[tmp]
        return(x)
    }, genes=rownames(rlist))
    groups <- groups[sapply(groups, function(x) length(x$tfmode))>=minsize]
    if (is.null(nullpw)) {
        if (per>0) {
			nullpw <- list(eset=sapply(1:per, function(i, rlist) sample(rlist), rlist=rlist[, 1]))
			colnames(nullpw$eset) <- 1:per
			rownames(nullpw$eset) <- rownames(rlist)
		}
		else nullpw <- ppwea3NULLf(groups)
    }
    if (names(nullpw)[1]=="eset") nullpw <- pwea3NULLgroups(nullpw, groups, cores=cores, verbose=verbose)
    if (names(nullpw)[1]=="groups") nullpw <- pwea3NULLf(nullpw, cores=cores, verbose=verbose)
    pb <- NULL
    if (verbose) {
        message("\nComputing enrichment by PWEA3 on ", length(groups), " gene-sets.")
        message("Process started at ", date())
    }
    if (cores>1) {
        es <- .viperParallelApply(1:length(groups), function(i, reg, r1, r2) {
            reg <- reg[[i]]
            pos <- match(names(reg$tfmode), rownames(r1))
            ll <- reg$likelihood / sum(reg$likelihood)
            sum1 <- t(filterRowMatrix(r2, pos)) %*% matrix(reg$tfmode * ll, length(ll), 1)
            ss <- sign(sum1)
            ss[ss==0] <- 1
            sum2 <- t(filterRowMatrix(r1, pos)) %*% matrix(ll * (1-abs(reg$tfmode)), length(ll), 1)
            return((abs(sum1)+(sum2*(sum2>0))) * ss)
        }, reg=groups, r1=rlist1, r2=rlist2, mc.cores=cores)
        es <- sapply(es, function(x) x)    
    }
    else {
        if (verbose) pb <- txtProgressBar(max=length(groups), style=3)
        es <- sapply(1:length(groups), function(i, reg, r1, r2, pb, verbose) {
            reg <- reg[[i]]
            pos <- match(names(reg$tfmode), rownames(r1))
            ll <- reg$likelihood / sum(reg$likelihood)
            sum1 <- t(filterRowMatrix(r2, pos)) %*% matrix(reg$tfmode * ll, length(ll), 1)
            ss <- sign(sum1)
            ss[ss==0] <- 1
            if (verbose) setTxtProgressBar(pb, i)
            sum2 <- t(filterRowMatrix(r1, pos)) %*% matrix(ll * (1-abs(reg$tfmode)), length(ll), 1)
            return((abs(sum1)+(sum2*(sum2>0))) * ss)
        }, reg=groups, r1=rlist1, r2=rlist2, pb=pb, verbose=verbose)
    }
    names(es) <- names(groups)
    if (is.vector(es)) es <- matrix(es, 1, length(es), dimnames=list(NULL, names(es)))
	colnames(es) <- names(groups)
    temp <- lapply(colnames(es), function(x, es, pfun, alter) pfun[[x]](es[, x]), es=es, pfun=nullpw)
    names(temp) <- colnames(es)
    nes <- sapply(temp, function(x) qnorm(x$p.value/2, lower.tail=FALSE))*sign(es)
	senes <- sqrt(fcvarna(nes)/nrow(nes))[, 1]
	if (nrow(nes)==1) senes <- NULL
	nes1 <- colMeans(nes)
    pval1 <- pnorm(abs(nes1), lower.tail=FALSE)*2
    if (verbose) message("\nProcess ended at ", date())
    return(list(es=colMeans(es), nes=nes1, nes.se=senes, size=sapply(groups, function(x) length(x$tfmode)), p.value=pval1, nes.bt=t(nes)))
}

#' Regulon-specific NULL model
#'
#' This function generates the regulon-specific NULL models
#'
#' @param pwnull Numerical matrix representing the null model, with genes as rows (geneID as rownames) and permutations as columns
#' @param groups List containing the regulons
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A list containing two elements:
#' \describe{
#' \item{groups}{Regulon-specific NULL model containing the enrichment scores}
#' \item{ss}{Direction of the regulon-specific NULL model}
#' }

pwea3NULLgroups <- function(pwnull, groups, cores=1, verbose=TRUE) {
	if (is.null(pwnull)) return(NULL)
    if (is.matrix(pwnull)) pwnull <- list(eset=pwnull)
    groups <- updateRegulon(groups)
    if (is.list(pwnull$eset)) {
        t <- unlist(pwnull$eset, use.names=FALSE)
        dim(t) <- c(length(pwnull$eset[[1]]), length(pwnull$eset))
        rownames(t) <- names(pwnull$eset[[1]])
        colnames(t) <- 1:length(pwnull$eset)
    }
    else t <- pwnull$eset
	t2 <- apply(t, 2, rank)/(nrow(t)+1)*2-1
	t1 <- abs(t2)*2-1
    t1 <- t1+(1-max(t1))/2
    t1 <- qnorm(t1/2+.5)
    t2 <- qnorm(t2/2+.5)
    pb <- NULL
    if (verbose) {
        message("\nComputing the null distribution for ", length(groups), " gene-sets.")
        message("Process started at ", date())
    }
    if (cores>1) {
        temp <- .viperParallelApply(1:length(groups), function(i, groups, t1, t2) {
            x <- groups[[i]]
            pos <- match(names(x$tfmode), rownames(t1))
            sum1 <- matrix(x$tfmode * x$likelihood, 1, length(x$tfmode)) %*% filterRowMatrix(t2, pos)
            ss <- sign(sum1)
            ss[ss==0] <- 1
            sum2 <- matrix((1-abs(x$tfmode)) * x$likelihood, 1, length(x$tfmode)) %*% filterRowMatrix(t1, pos)
            return(list(es=as.vector(abs(sum1) + sum2*(sum2>0)) / sum(x$likelihood), ss=ss))
        }, groups=groups, t1=t1, t2=t2, mc.cores=cores)
    }
    else {
        if (verbose) pb <- txtProgressBar(max=length(groups), style=3)
        temp <- lapply(1:length(groups), function(i, groups, t1, t2, pb, verbose) {
            x <- groups[[i]]
            pos <- match(names(x$tfmode), rownames(t1))
            sum1 <- matrix(x$tfmode * x$likelihood, 1, length(x$tfmode)) %*% filterRowMatrix(t2, pos)
            ss <- sign(sum1)
            ss[ss==0] <- 1
            if (verbose) setTxtProgressBar(pb, i)
            sum2 <- matrix((1-abs(x$tfmode)) * x$likelihood, 1, length(x$tfmode)) %*% filterRowMatrix(t1, pos)
            return(list(es=as.vector(abs(sum1) + sum2*(sum2>0)) / sum(x$likelihood), ss=ss))
        }, groups=groups, pb=pb, t1=t1, t2=t2, verbose=verbose)
    }
    names(temp) <- names(groups)
    if (verbose) message("\nProcess ended at ", date())
    es <- t(sapply(temp, function(x) x$es))
    ss <- t(sapply(temp, function(x) x$ss))
    return(list(groups=es, ss=ss))
}

#' Null model function
#'
#' This function generates the NULL model function, which computes the normalized enrichment score and associated p-value
#'
#' @param pwnull Object generated by \code{pwea3NULLgroups} function
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return List of function to compute NES and p-value

pwea3NULLf <- function(pwnull, cores=1, verbose=TRUE) {
	if (!is.null(pwnull)) {
        pb <- NULL
        if (verbose) {
            message("\nComputing the null distribution function for ", nrow(pwnull$groups), " TFs.")
            message("Process started at ", date())
        }
        if (cores>1) {
            res <- .viperParallelApply(1:nrow(pwnull$groups), function(i, pwnull) {
                x <- pwnull$groups[i, ]
                iqr <- quantile(x, c(.5, 1-5/length(x)))
                pd <- ecdf(x)
                a <- list(x=knots(pd), y=pd(knots(pd)))
                filtro <- a$x>iqr[1] & a$x<=iqr[2] & a$y<1
                spl <- smooth.spline(a$x[filtro], -log(1-a$y[filtro]), spar=.75)
                return(function(x, alternative="greater") {
                    x <- abs(x)
                    p <- exp(-predict(spl, x)$y)
                    pos <- which(x<iqr[1])
                    if (length(pos)>0) p[pos] <- 1-pd(x[pos])
                    nes <- qnorm(p/2, lower.tail=FALSE)
                    if (length(p)==1) {
                        p <- as.vector(p)
                        nes <- as.vector(nes)
                    }
                    list(nes=nes*pwnull$ss[i], p.value=p)
                })
            }, pwnull=pwnull, mc.cores=cores)
        }
        else {
            if (verbose) pb <- txtProgressBar(max=nrow(pwnull$groups), style=3)
            res <- lapply(1:nrow(pwnull$groups), function(i, pwnull, pb, verbose) {
                if (verbose) setTxtProgressBar(pb, i)
                x <- pwnull$groups[i, ]
                iqr <- quantile(x, c(.5, 1-5/length(x)))
                pd <- ecdf(x)
                a <- list(x=knots(pd), y=pd(knots(pd)))
                filtro <- a$x>iqr[1] & a$x<=iqr[2] & a$y<1
                spl <- smooth.spline(a$x[filtro], -log(1-a$y[filtro]), spar=.75)
                return(function(x, alternative="greater") {
                    x <- abs(x)
                    p <- exp(-predict(spl, x)$y)
                    pos <- which(x<iqr[1])
                    if (length(pos)>0) p[pos] <- 1-pd(x[pos])
                    nes <- qnorm(p/2, lower.tail=FALSE)
                    if (length(p)==1) {
                        p <- as.vector(p)
                        nes <- as.vector(nes)
                    }
                    list(nes=nes*pwnull$ss[i], p.value=p)
                })
            }, pwnull=pwnull, pb=pb, verbose=verbose)
        }
        if (verbose) message("\nProcess ended at ", date(), "\n", sep="")
		names(res) <- rownames(pwnull$groups)
		return(res)
	}
}

ppwea3NULLf <- function(regulon) {
    lapply(regulon, function(x) {
        ww <- x$likelihood
        ww <- ww/max(ww)
        ww <- sqrt(sum(ww^2))
        return(function(x, alternative="two.sided") {
            x <- x*ww
            p <- switch(pmatch(alternative, c("two.sided", "less", "greater")),
                pnorm(abs(x), lower.tail=FALSE)*2,
                pnorm(x, lower.tail=TRUE),
                pnorm(x, lower.tail=FALSE))
        list(nes=x, p.value=p)
        })
    })
}

regulonScore <- function(x, slope=20, inflection=.3) {
    if (is.list(x)) return(lapply(x, regulonScore, slope=slope, inflection=inflection))
    return(1/(1+inflection^(slope*(abs(x)-inflection)))*sign(x))
}

regulonDist <- function (groups, cutoff = 50, method = "FET", copula = TRUE, mode = FALSE) {
    if (!mode) 
        groups <- lapply(groups, function(x) {
            res <- abs(sign(x))
            res[res == 0] <- 1
            return(res)
        })
    total <- unique(unlist(lapply(groups, names), use.names = FALSE))
    groups <- groups[sapply(groups, length) >= cutoff]
    group <- groups
    dmatrix <- NULL
    for (i in 1:(length(groups) - 1)) {
        group <- group[-which(names(group) == names(groups)[i])]
        test <- groups[[i]]
        res <- sapply(group, function(x, test, total, method) {
            tmp1 <- tmp2 <- rep(0, length(total))
            tmp1[match(names(test), total)] <- sign(test)
            tmp2[match(names(x), total)] <- sign(x)
            tmp <- abs(tmp1 + tmp2) > 0
            test1 <- abs(tmp1 * tmp) > 0
            test2 <- abs(tmp2 * tmp) > 0
            tmp <- fisher.test(table(test1, test2), alternative = "greater")
            switch(pmatch(method, c("odds", "FET")), res <- tmp$estimate, 
                res <- tmp$p.value)
            res
        }, test = test, total = total, method = method)
        dmatrix <- cbind(dmatrix, c(dmatrix[i, ], 0, res))
    }
    dmatrix <- cbind(dmatrix, c(dmatrix[length(groups), ], 0))
    rownames(dmatrix) <- colnames(dmatrix) <- names(groups)
    if (copula) {
        rdm <- rank(dmatrix)
        dmatrix <- matrix(rdm/max(rdm), nrow(dmatrix), ncol(dmatrix), 
            dimnames = list(names(groups), names(groups)))
    }
    return(dmatrix)
}

#' Filter for rows of a matrix with no loss of col and row names
#'
#' This function filters the rows of a matrix returning always a two dimensional matrix
#'
#' @param x Matrix
#' @param filter Logical or numerical index of rows
#' @return Matrix
#' @export
filterRowMatrix <- function(x, filter) {
    if (is.logical(filter)) largo <- length(which(filter))
    else largo <- length(filter)
    matrix(x[filter, ], largo, ncol(x), dimnames=list(rownames(x)[filter], colnames(x)))
}

#' Filter for columns of a matrix with no loss of col and row names
#'
#' This function filters the columns of a matrix returning always a two dimensional matrix
#'
#' @param x Matrix
#' @param filter Logical or numerical index of columns
#' @return Matrix
#' @export
filterColMatrix <- function(x, filter) t(filterRowMatrix(t(x), filter))

#' Mode of continuous distributions
#' 
#' This function computes the mode for continuous distributions
#' 
#' @param x Numeric data vector
#' @param adj Number indicating the adjustment for the kernel bandwidth
#' @return Number
#' @examples
#' data(bcellViper, package="bcellViper")
#' d1 <- exprs(dset)
#' mean(d1[, 1])
#' median(d1[, 1])
#' distMode(d1[, 1])
#' plot(density(d1[, 1]))
#' abline(v=c(mean(d1[, 1]), median(d1[, 1]), distMode(d1[, 1])), col=c("green", "red", "blue"))
#' legend("topleft", c("Mean", "Median", "Mode"), col=c("green", "red", "blue"), lwd=4)
#' @export
distMode <- function (x, adj = 1) {
    tmp <- density(x, adjust = adj)
    return(tmp$x[which.max(tmp$y)])
}

#' Correction for pleiotropy
#' 
#' This function penalyze the regulatory interactions based on pleiotropy analysis
#' 
#' @param ss Named vector containing the gene expression signature
#' @param nes Named vector containing the normalized enrichment scores
#' @param regul Regulon object
#' @param regulators Number indicating the number of top regulators to consider for the analysis or the p-value threshold for considering significant regulators
#' @param shadow Number indicating the p-value threshold for considering a significant shadow effect
#' @param targets Integer indicating the minimal number of overlaping targets to consider a pair of regulators for pleiotropy analysis
#' @param penalty Number higher than 1 indicating the penalty for the pleiotropic interactions. 1 = no penalty
#' @param method Character string indicating the method to use for computing the pleiotropy, either absolute or adaptive
#' @return Corrected regulon object

shadowRegulon <- function(ss, nes, regul, regulators=.05, shadow=.05, targets=10, penalty=2, method=c("absolute", "adaptive")) {
    method <- match.arg(method)
    pval <- pnorm(abs(nes), lower.tail=FALSE)*2
    if (regulators<1) tfs <- names(pval)[pval<regulators]
    else tfs <- names(pval)[order(pval)[1:regulators]]
    pos <- grep("--", tfs)
    if (length(pos)>0) tfs <- tfs[-pos]
    if (length(tfs)<2) return(NULL)
    # Cache each directed target overlap once. The same intersections are used
    # to construct shadow regulons and later to apply their penalties.
    tfTargets <- lapply(tfs, function(tf) names(regul[[tf]]$tfmode))
    names(tfTargets) <- tfs
    overlaps <- lapply(tfs, function(tf1) {
        others <- tfs[tfs != tf1]
        out <- lapply(others, function(tf2) which(tfTargets[[tf1]] %in% tfTargets[[tf2]]))
        names(out) <- others
        out
    })
    names(overlaps) <- tfs
    tmp <- lapply(unique(tfs), function(tf1, tfs, regul, ss, nes, targets) {
        reg <- lapply(tfs[tfs != tf1], function(tf2, regul, tf1) {
            pos <- overlaps[[tf1]][[tf2]]
            list(tfmode=regul[[tf1]]$tfmode[pos], likelihood=regul[[tf1]]$likelihood[pos])
        }, regul=regul, tf1=tf1)
        names(reg) <- tfs[tfs != tf1]
        pos <- which(names(ss) %in% names(regul[[tf1]]$tfmode))
        s2 <- rank(ss[pos])/(length(ss[pos])+1)*2-1
        s1 <- abs(s2)*2-1
        s1 <- s1+(1-max(s1))/2
        s1 <- qnorm(s1/2+.5)
        tmp <- sign(nes[tf1])
        if (tmp==0) tmp <- 1
        s2 <- qnorm(s2/2+.5)*tmp
        tmp <- sapply(reg, function(x, s1, s2, targets) {
            if (length(x$tfmode)<targets) return(NA)
            pos <- match(names(x$tfmode), names(s1))
            sum1 <- sum(x$tfmode * x$likelihood * s2[pos])
            ss <- sign(sum1)
            ss[ss==0] <- 1
            sum2 <- sum((1-abs(x$tfmode)) * x$likelihood * s1[pos])
            ww <- x$likelihood/max(x$likelihood)
            return((abs(sum1) + sum2*(sum2>0)) / sum(x$likelihood) * sign(ss) * sqrt(sum(ww^2)))
        }, s1=s1, s2=s2, targets=targets)
        return(pnorm(tmp, lower.tail=FALSE))
    }, tfs=tfs, regul=regul, ss=ss, nes=nes, targets=targets)
    names(tmp) <- unique(tfs)
    pval <- unlist(tmp, use.names=F)
    names(pval) <- paste(rep(names(tmp), sapply(tmp, length)), unlist(lapply(tmp, names), use.names=FALSE), sep=" x ")
    pval <- pval[!is.na(pval)]
    regind <- t(combn(tfs, 2))
    regind <- filterRowMatrix(regind, paste(regind[, 1], regind[, 2], sep=" x ") %in% names(pval))
    pval <- cbind(pval[match(paste(regind[, 1], regind[, 2], sep=" x "), names(pval))], pval[match(paste(regind[, 2], regind[, 1], sep=" x "), names(pval))])
    tests <- table(as.vector(regind))
    switch(method,
    absolute={
        tmp <- rbind(regind[pval[, 1]<shadow & pval[, 2]>shadow, ], regind[, 2:1][pval[, 1]>shadow & pval[, 2]<shadow, ])
        if (nrow(tmp)==0) return(NULL)
        for (i in 1:nrow(tmp)) {
            ll <- regul[[tmp[i, 1]]]$likelihood
            pos <- overlaps[[tmp[i, 1]]][[tmp[i, 2]]]
            ll[pos] <- ll[pos]/penalty^(1/tests[tmp[i, 1]])
            regul[[tmp[i, 1]]]$likelihood <- ll
        }
    },
    adaptive={
        pval1 <- log10(pval[, 2])-log10(pval[, 1])
        tmp <- NULL
        if (length(which(pval1>0))>0) tmp <- filterRowMatrix(regind, pval1>0)
        if (length(which(pval1<0))>0) tmp <- rbind(tmp, filterRowMatrix(regind, pval1<0)[, 2:1])
        pval1 <- c(pval1[pval1>0], -pval1[pval1<0])
        if (is.null(nrow(tmp))) return(NULL)
        for (i in 1:nrow(tmp)) {
            ll <- regul[[tmp[i, 1]]]$likelihood
            pos <- overlaps[[tmp[i, 1]]][[tmp[i, 2]]]
            ll[pos] <- ll[pos]/(1+pval1[i])^(penalty/tests[tmp[i, 1]])
            regul[[tmp[i, 1]]]$likelihood <- ll
        }
    })
    return(regul[which(names(regul) %in% tmp[, 1])])
}

#' Approximate empirical commulative distribution function
#'
#' This function generates an empirical null model that computes a normalized statistics and p-value
#' 
#' @param dnull Numerical vector representing the null model
#' @param symmetric Logical, whether the distribution should betreated as symmetric around zero and only one tail should be approximated
#' @param n Integer indicating the number of points to evaluate the empirical cummulative probability function
#' @return function with two parameters, \code{x} and \code{alternative}

aecdf <- function(dnull, symmetric=FALSE, n=100) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        tmp <- sort(abs(dnull), decreasing=T)
        i <- 4
        n <- 4
        while(n<14) {
            i <- i+1
            n <- length(unique(tmp[1:i]))
            if (n==5) iq1 <- i
        }
        tl1 <- i
        iqr <- quantile(abs(dnull), c(.5, 1-iq1/length(dnull)))
        epd <- ecdf(abs(dnull))
        a <- list(x=abs(dnull), y=epd(abs(dnull)))
        a1 <- list(x=a$x[length(a$x)-(tl1:iq1)+1]-iqr[2], y=log(1-epd(iqr[2]))-log(1-a$y[length(a$x)-(tl1:iq1)+1]))
        a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
        if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
        fit <- lm(y~0+x, data=a1)
        val <- seq(0, iqr[2], length=n)
        pd <- approxfun(val, epd(val), method="linear", yleft=0, rule=2)
        dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
            alternative <- match.arg(alternative)
            x1 <- abs(x)
            p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
            p[!is.finite(p)] <- 1
            p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
            nes <- qnorm(p/2, lower.tail=F)*sign(x)
            switch(alternative,
                   two.sided={p <- p},
                   greater={p <- p/2; p[x<0] <- 1-p[x<0]},
                   less={p <- p/2; p[x>0] <- 1-p[x>0]}
            )
            names(nes) <- names(p) <- names(x)
            list(nes=nes, p.value=p)
        }
        return(dnull)
    }
    tmp <- sort(dnull, decreasing=FALSE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq1 <- i
    }
    tl1 <- i
    tmp <- sort(dnull, decreasing=TRUE)
    i <- 4
    n <- 4
    while(n<14) {
        i <- i+1
        n <- length(unique(tmp[1:i]))
        if (n==5) iq2 <- i
    }
    tl2 <- i
    iqr <- quantile(dnull, c(iq1/length(dnull), .5, 1-iq2/length(dnull)))
    epd <- ecdf(dnull)
    a <- list(x=dnull, y=epd(dnull))
    a1 <- list(x=a$x[iq1:tl1]-iqr[1], y=log(epd(iqr[1]))-log(a$y[iq1:tl1]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
    if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
    fit1 <- lm(y~0+x, data=a1)
    a1 <- list(x=a$x[length(a$x)-(tl2:iq2)+1]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(tl2:iq2)+1]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
    if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
    fit2 <- lm(y~0+x, data=a1)
    val <- seq(iqr[1], iqr[3], length=n)
    pd <- approxfun(val, epd(val), method="linear", rule=2)
    dnull <- function(x, alternative=c("two.sided", "greater", "less")) {
        alternative <- match.arg(alternative)
        p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
        p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
        p1[!is.finite(p1)] <- 1
        p2[!is.finite(p2)] <- 1
        p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
        nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
        switch(alternative,
               two.sided={p <- p*2},
               greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
               less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
        )
        names(nes) <- names(p) <- names(x)
        list(nes=nes, p.value=p)
    }
    return(dnull)
}

aecdf1 <- function(dnull, symmetric=FALSE, x, alternative=c("two.sided", "greater", "less")) {
    dnull <- dnull[is.finite(dnull)]
    if (symmetric) {
        iqr <- quantile(abs(dnull), c(.5, 1-5/length(dnull)))
        pd <- ecdf(abs(dnull))
        a <- list(x=abs(dnull), y=pd(abs(dnull)))
        a1 <- list(x=a$x[length(a$x)-(tl2:iq2)+1]-iqr[3], y=log(1-epd(iqr[3]))-log(1-a$y[length(a$x)-(tl2:iq2)+1]))
        a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
        if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
        a1 <- list(x=a$x[length(a$x)-(15:4)]-iqr[2], y=log(1-pd(iqr[2]))-log(1-a$y[length(a$x)-(15:4)]))
        a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
        if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
        fit <- lm(y~0+x, data=a1)
        alternative <- match.arg(alternative)
        x1 <- abs(x)
        p <- exp(log(1-pd(iqr[2]))-predict(fit, list(x=x1-iqr[2])))
        p <- p * (x1>iqr[2]) + (1-pd(x1)) * (x1<=iqr[2])
        nes <- qnorm(p/2, lower.tail=FALSE)*sign(x)
        switch(alternative,
               two.sided={p <- p},
               greater={p <- p/2; p[x<0] <- 1-p[x<0]},
               less={p <- p/2; p[x>0] <- 1-p[x>0]}
        )
        names(nes) <- names(p) <- names(x)
        return(list(nes=nes, p.value=p))
    }
    iqr <- quantile(dnull, c(5/length(dnull), .5, 1-5/length(dnull)))
    pd <- ecdf(dnull)
    a <- list(x=dnull, y=pd(dnull))
    a1 <- list(x=a$x[5:14]-iqr[1], y=log(pd(iqr[1]))-log(a$y[5:14]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
    if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
    fit1 <- lm(y~0+x, data=a1)
    a1 <- list(x=a$x[length(a$x)-(15:4)]-iqr[3], y=log(1-pd(iqr[3]))-log(1-a$y[length(a$x)-(15:4)]))
    a1 <- lapply(a1, function(x, pos) x[pos], pos=which(is.finite(a1$y)))
    if (length(a1$x)<3) stop("Not enough permutations to compute NULL distribution", call.=FALSE)
    fit2 <- lm(y~0+x, data=a1)
    alternative <- match.arg(alternative)
    p1 <- exp(log(pd(iqr[1]))-predict(fit1, list(x=x-iqr[1])))
    p2 <- exp(log(1-pd(iqr[3]))-predict(fit2, list(x=x-iqr[3])))
    p <- p1*(x<iqr[1]) + p2*(x>iqr[3]) + pd(x)*(x>=iqr[1] & x<iqr[2]) + (1-pd(x))*(x>=iqr[2] & x<=iqr[3])
    nes <- qnorm(p, lower.tail=F)*sign(x-iqr[2])
    switch(alternative,
           two.sided={p <- p*2},
           greater={p[x<iqr[2]] <- 1-p[x<iqr[2]]},
           less={p[x>=iqr[2]] <- 1-p[x>=iqr[2]]}
    )
    names(nes) <- names(p) <- names(x)
    return(list(nes=nes, p.value=p))
}

#' Sigmoid transformation
#' 
#' This function transforms a numeric vector using a sigmoid function
#' 
#' @param x Numeric vector
#' @param slope Number indicating the slope at the inflection point
#' @param inflection Number indicating the inflection point
#' @return Numeric vector
sigT <- function (x, slope = 20, inflection = 0.5) 1 - 1/(1 + exp(slope * (x - inflection)))
