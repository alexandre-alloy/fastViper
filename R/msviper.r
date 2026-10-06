#' msVIPER
#'
#' This function performs MAster Regulator INference Analysis
#'
#' @param ges Vector containing the gene expression signature to analyze, or matrix with columns containing bootstraped signatures
#' @param regulon Object of class regulon
#' @param nullmodel Matrix of genes by permutations containing the NULL model signatures. A parametric approach equivalent to shuffle genes will be used if nullmodel is ommitted.
#' @param pleiotropy Logical, whether correction for pleiotropic regulation should be performed
#' @param minsize Number indicating the minimum allowed size for the regulons
#' @param adaptive.size Logical, whether the weight (likelihood) should be used for computing the regulon size
#' @param ges.filter Logical, whether the gene expression signature should be limited to the genes represented in the interactome
#' @param synergy Number indicating the synergy computation mode: (0) for no synergy computation; (0-1) for establishing the p-value cutoff for individual TFs to be included in the synergy analysis; (>1) number of top TFs to be included in the synergy analysis
#' @param level Integer, maximum level of combinatorial regulation
#' @param pleiotropyArgs list of 5 numbers for the pleotropy correction indicating: regulators p-value threshold, pleiotropic interaction p-value threshold, minimum number of targets in the overlap between pleiotropic regulators, penalty for the pleiotropic interactions and the pleiotropy analysis method, either absolute or adaptive
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A msviper object containing the following components:
#' \describe{
#' \item{signature}{The gene expression signature}
#' \item{regulon}{The final regulon object used}
#' \item{es}{Enrichment analysis results including regulon size, normalized enrichment score and p-value}
#' \item{param}{msviper parameters, including \code{minsize}, \code{adaptive.size}}
#' }
#' @seealso \code{\link{viper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- rowTtest(dset, "description", c("CB", "CC"), "N")$statistic
#' dnull <- ttestNull(dset, "description", c("CB", "CC"), "N", per=100) # Only 100 permutations to reduce computation time, but it is recommended to perform at least 1000 permutations
#' mra <- msviper(sig, regulon, dnull)
#' plot(mra, cex=.7)
#' @export

msviper <- function(ges, regulon, nullmodel=NULL, pleiotropy=FALSE, minsize=25, adaptive.size=FALSE, ges.filter=TRUE, synergy=0, level=10, pleiotropyArgs=list(regulators=.05, shadow=.05, targets=10, penalty=20, method="adaptive"), cores=1, verbose=TRUE) {
# Cleaning the parameters
	if (is.vector(ges)) ges <- matrix(ges, length(ges), 1, dimnames=list(names(ges), NULL))
	regulon <- updateRegulon(regulon)
	if (ges.filter) ges <- filterRowMatrix(ges, rownames(ges) %in% unique(c(names(regulon), unlist(lapply(regulon, function(x) names(x$tfmode)), use.names=FALSE))))
	regulon <- lapply(regulon, function(x, genes) {
        filtro <- names(x$tfmode) %in% genes
        x$tfmode <- x$tfmode[filtro]
        x$likelihood <- x$likelihood[filtro]
        return(x)
    }, genes=rownames(ges))
	if (!is.null(nullmodel)) nullmodel <- nullmodel[match(rownames(ges), rownames(nullmodel)), ]
# Running msviper
	if (adaptive.size) {
        regulon <- regulon[sapply(regulon, function(x) {
            sum(x$likelihood/max(x$likelihood))
        })>=minsize]
	}
	else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode))>=minsize]
    if (verbose) message("Computing regulon enrichment with aREA algorithm")
    res <- aREA(ges, regulon, method="loop", cores=cores, minsize=0, verbose=verbose)
    if (is.null(nullmodel)) nes <- res$nes
    else {
        if (verbose) message("\nEstimating the normalized enrichment scores")
        tmp <- aREA(nullmodel, regulon, cores=cores, minsize=0, verbose=verbose)$es
        nes <- t(sapply(1:nrow(tmp), function(i, tmp, es) {
            aecdf(tmp[i, ], symmetric=TRUE)(es[i, ])$nes
        }, tmp=tmp, es=res$es))
        if (nrow(nes)==1) nes <- t(nes)
        rownames(nes) <- rownames(res$es)
    }
    if (pleiotropy) {
        pb <- NULL
        if (verbose) {
            message("\nComputing pleiotropy for ", ncol(nes), " samples.")
            message("\nProcess started at ", date())
        }
        if (cores>1) {
            nes <- .viperParallelApply(1:ncol(nes), function(i, ss, nes, regulon, args, dnull) {
                nes <- nes[, i]
                sreg <- shadowRegulon(ss[, i], nes, regulon, regulators=args[[1]], shadow=args[[2]], targets=args[[3]], penalty=args[[4]], method=args[[5]])
                if (!is.null(sreg)) {
                    if (is.null(dnull)) tmp <- aREA(ss[, i], sreg, minsize=5, cores=1)$nes[, 1]
                    else {
                        tmp <- aREA(cbind(ss[, i], dnull), sreg, minsize=5, cores=1)$es
                        tmp <- apply(tmp, 1, function(x) aecdf(x[-1], symmetric=TRUE)(x[1])$nes)
                    }
                    nes[match(names(tmp), names(nes))] <- tmp
                }
                return(nes)
            }, ss=ges, nes=nes, regulon=regulon, args=pleiotropyArgs, dnull=nullmodel, mc.cores=cores)
            nes <- sapply(nes, function(x) x)    
        }
        else {
            if (verbose) pb <- txtProgressBar(max=ncol(nes), style=3)
            nes <- sapply(1:ncol(nes), function(i, ss, nes, regulon, args, dnull, pb) {
                nes <- nes[, i]
                sreg <- shadowRegulon(ss[, i], nes, regulon, regulators=args[[1]], shadow=args[[2]], targets=args[[3]], penalty=args[[4]], method=args[[5]])
                if (!is.null(sreg)) {
                    if (is.null(dnull)) tmp <-aREA(ss[, i], sreg, minsize=5)$nes[, 1]
                    else {
                        tmp <- aREA(cbind(ss[, i], dnull, minsize=5), sreg)$es
                        tmp <- apply(tmp, 1, function(x) aecdf(x[-1], symmetric=TRUE)(x[1])$nes)
                    }
                    nes[match(names(tmp), names(nes))] <- tmp
                }
                if (is(pb, "txtProgressBar")) setTxtProgressBar(pb, i)    
                return(nes)
            }, ss=ges, nes=nes, regulon=regulon, args=pleiotropyArgs, dnull=nullmodel, pb=pb)
        }
        if (verbose) message("\nProcess ended at ", date(), "\n")
        colnames(nes) <- colnames(ges)
    }
    senes <- sqrt(frvarna(nes)/ncol(nes))[, 1]
    if (ncol(nes)==1) senes <- NULL
    nes1 <- rowMeans(nes)
    pval1 <- pnorm(abs(nes1), lower.tail=FALSE)*2
    res <- list(signature=ges, regulon=regulon, es=list(nes=nes1, nes.se=senes, size=sapply(regulon, function(x) length(x$tfmode)), p.value=pval1, nes.bt=nes), param=list(minsize=minsize, adaptive.size=adaptive.size), nullmodel=nullmodel)
    class(res) <- "msviper"
	if (synergy>0) res <- msviperSynergy(msviperCombinatorial(mobj=res, regulators=synergy, nullmodel=nullmodel, minsize=minsize, adaptive.size=adaptive.size, level=level, verbose=verbose), verbose=verbose)
	return(res)
}

#' msviper bootstraps integration
#' 
#' This function integrates the bootstrap msviper results
#' 
#' @param mobj msviper object
#' @param method Character string indicating the method to use, either mean, median or mode
#' @return msviper object
#' @seealso \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- bootstrapTtest(dset, "description", c("CB", "CC"), "N")
#' mra <- msviper(sig, regulon)
#' plot(mra, cex=.7)
#' @export
bootstrapmsviper <- function(mobj, method=c("mean", "median", "mode")) {
    method <- match.arg(method)
    if (ncol(mobj$es$nes.bt)>1) switch(method,
        mean={mobj$es$nes <- rowMeans(mobj$es$nes.bt)},
        median={mobj$es$nes <- apply(mobj$es$nes.bt, 1, median)},
        mode={mobj$es$nes <- apply(mobj$es$nes.bt, 1, distMode)})
    mobj$es$p.value <- pnorm(abs(mobj$es$nes), lower.tail=FALSE)*2
    return(mobj)
}

#' msviper combinatorial analysis
#'
#' This function performs combinatorial analysis for msviper objects
#'
#' @param mobj msviper object generated by \code{msviper} function
#' @param regulators Either a number between 0 and 1 indicating the p-value cutoff for individual TFs to be included in the combinations analysis; (>1) indicating the number of top TFs to be included in the combinations analysis; or a vector of character strings indicating the TF IDs to be included in the analysis
#' @param nullmodel Matrix of genes by permutations containing the NULL model signatures. Taken from \code{mobj} by default
#' @param minsize Number indicating the minimum allowed size for the regulons, taken from \code{mobj} by default
#' @param adaptive.size Logical, whether the weight (likelihood) should be used for computing the size, taken from \code{mobj} by default
#' @param level Integer, maximum level of combinatorial regulation
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param processAll Logical, whether all pairs, even if not significant, should be processed for synergy
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return A msviper object
#' @seealso \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- rowTtest(dset, "description", c("CB", "CC"), "N")$statistic
#' dnull <- ttestNull(dset, "description", c("CB", "CC"), "N", per=100) # Only 100 permutations to reduce computation time, but it is recommended to perform at least 1000 permutations
#' mra <- msviper(sig, regulon, dnull)
#' mra <- msviperCombinatorial(mra, 20)
#' plot(mra, cex=.7)
#' @export
msviperCombinatorial <- function(mobj, regulators=100, nullmodel=NULL, minsize=NULL, adaptive.size=NULL, level=10, cores=1, processAll=FALSE, verbose=TRUE) {
    synergy <- regulators
    if (is.null(minsize)) minsize <- mobj$param$minsize
	if (is.null(adaptive.size)) adaptive.size <- mobj$param$adaptive.size
    if (is.null(nullmodel)) nullmodel <- mobj$nullmodel
    if (length(synergy)==1 & synergy[1]==0) return(mobj)
	res <- mobj$es
	regulon <- mobj$regulon
    tfs <- synergy
    if (length(synergy)==1) {
        if (synergy<1) tfs <- names(res$p.value)[res$p.value<synergy]
        else tfs <- names(res$p.value)[order(res$p.value)[1:synergy]]
    }
    res1 <- list(es=res)
	resp <- res$nes
	res <- list(res)
	regul <- NULL
	n=2
	while(length(tfs)>n & n<=level) {
		if (verbose) message("\n-------------------------------------------\nComputing synergy for combination of ", n, " TFs\n-------------------------------------------\n\n")
		res1 <- comregulationAnalysis(combn(tfs, n), mobj$signature, regulon, nullmodel, resp, res1$es$p.value, minsize, adaptive.size, cores=cores, processAll=processAll, verbose=verbose)
		if(length(res1)==0) break
		res <- c(res, list(res1$es))
		regul <- c(regul, res1$regul)
		n <- n+1
		tfs <- unique(unlist(strsplit(as.character(names(res1$es$nes)), "--"), use.names=FALSE))
	}	
	names(res) <- 1:length(res)
	res <- res[sapply(res, function(x) length(x$nes))>0]
	tmp <- names(res[[1]])
	res <- lapply(tmp, function(nom, res) {
		tmp <- unlist(lapply(res, function(x, nom) x[[nom]], nom=nom), use.names=FALSE)
		names(tmp) <- unlist(lapply(res, function(x, nom) names(x[[nom]]), nom=nom), use.names=FALSE)
		return(tmp)
	}, res=res)
	names(res) <- tmp
	regulon <- c(regulon, regul)
	mobj$regulon <- regulon
	mobj$es <- res
    class(mobj) <- "msviper"
	return(mobj)
}

#' msviper synergy analysis
#'
#' This function performs a synergy analysis for combinatorial regulation
#'
#' @param mobj msviper object containing combinatorial regulation results generated by \code{msviperCombinatorial}
#' @param per Integer indicating the number of permutations
#' @param seed Integer indicating the seed for the permutations, 0 for disable it
#' @param cores Integer number of worker processes. Values above 1 use PSOCK on Windows and fork workers on Unix.
#' @param verbose Logical, whether progression messages should be printed in the terminal
#' @return Updated msviper object containing the sygergy p-value
#' @seealso \code{\link{msviper}}
#' @examples
#' data(bcellViper, package="bcellViper")
#' sig <- rowTtest(dset, "description", c("CB", "CC"), "N")$statistic
#' dnull <- ttestNull(dset, "description", c("CB", "CC"), "N", per=100) # Only 100 permutations to reduce computation time, but it is recommended to perform at least 1000 permutations
#' mra <- msviper(sig, regulon, dnull)
#' mra <- msviperCombinatorial(mra, 20)
#' mra <- msviperSynergy(mra)
#' summary(mra)
#' @export

msviperSynergy <- function(mobj, per=1000, seed=1, cores=1, verbose=TRUE) {
    if (length(grep("--", names(mobj$regulon)))==0) stop("No significant co-regulons to evaluate", call.=FALSE)
    if (seed>0) set.seed(round(seed))
    pos <- which(sapply(strsplit(names(mobj$regulon), "--"), length)>1)
    pb <- NULL
    if (verbose) {
        message("\nComputing synergy statistics for ", length(pos), " co-regulons.")
        message("Process started at ", date())
    }
    if (cores>1) {
        res <- .viperParallelApply(1:length(pos), function(i, pos, ss, nes, regul) {
            i <- pos[i]
            pos <- match(unlist(strsplit(names(regul)[i], "--")), names(regul))
            tgs <- unique(unlist(lapply(regul[pos], function(x) names(x$tfmode)), use.names=FALSE))
            tfmode <- sapply(regul[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs)
            rownames(tfmode) <- tgs
            tfmode <- t(t(tfmode)*sign(nes[match(names(regul)[pos], names(nes))]))
            tfmode[is.na(tfmode)] <- 0
            ll <- sapply(regul[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs)
            ll[is.na(ll)] <- 0
            ll1 <- ll/rowSums(ll)
            tfmode <- rowSums(tfmode*ll1)
            ll <- apply(ll, 1, max)
            tgs <- names(regul[[i]]$tfmode)
            r2 <- apply(ss, 2, rank)/(nrow(ss)+1)*2-1
            r1 <- abs(r2)*2-1
            r1[r1==(-1)] <- 1-(1/length(r1))
            r1 <- qnorm(r1/2+.5)
            r2 <- qnorm(r2/2+.5)
            dim(r1) <- dim(r2) <- dim(ss)
            rownames(r1) <- rownames(r2) <- rownames(ss)
            pos <- match(names(tfmode), rownames(r1))
            r1 <- filterRowMatrix(r1, pos)
            r2 <- filterRowMatrix(r2, pos)
            dnull <- sapply(1:per, function(i, ll, samp) {
                ll[sample(length(ll), samp)] <- 0
                return(ll)
            }, ll=ll, samp=length(tgs))
            dnull <- cbind(ll, dnull)
            dnull[match(tgs, names(tfmode)), 1] <- 0
            dnull <- t(t(dnull)/colSums(dnull))		
            sum1 <- t(r2) %*% (dnull*tfmode)
            sum2 <- t(r1) %*% (dnull * (1-abs(tfmode)))
            es <- colMeans(abs(sum1)+sum2*(sum2>0))
            x <- es[1]
            es <- es[-1]
            iqr <- quantile(es, c(.5, 5/length(es)))
            pd <- ecdf(es)
            a <- list(x=knots(pd), y=pd(knots(pd)))
            filtro <- a$x<iqr[1] & a$x>=iqr[2] & a$y<1
            spl <- smooth.spline(a$x[filtro], -log(a$y[filtro]), spar=.75)
            p <- exp(-predict(spl, x)$y)
            pos <- which(x>iqr[1])
            if (x>iqr[1]) p <- pd(x)
            return(p)
        }, ss=mobj$signature, nes=mobj$es$nes, regul=mobj$regulon, pos=pos, mc.cores=cores)
        res <- sapply(res, function(x) x)    
    }
    else {
        if (verbose) pb <- txtProgressBar(max=length(pos), style=3)
        res <- sapply(1:length(pos), function(i, pos, ss, nes, regul, pb, verbose) {
            if (verbose) setTxtProgressBar(pb, i)
            i <- pos[i]
            pos <- match(unlist(strsplit(names(regul)[i], "--")), names(regul))
    		tgs <- unique(unlist(lapply(regul[pos], function(x) names(x$tfmode)), use.names=FALSE))
    		tfmode <- sapply(regul[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs)
    		rownames(tfmode) <- tgs
    		tfmode <- t(t(tfmode)*sign(nes[match(names(regul)[pos], names(nes))]))
    		tfmode[is.na(tfmode)] <- 0
    		ll <- sapply(regul[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs)
    		ll[is.na(ll)] <- 0
    		ll1 <- ll/rowSums(ll)
    		tfmode <- rowSums(tfmode*ll1)
    		ll <- apply(ll, 1, max)
    		tgs <- names(regul[[i]]$tfmode)
            r2 <- apply(ss, 2, rank)/(nrow(ss)+1)*2-1
            r1 <- abs(r2)*2-1
            r1[r1==(-1)] <- 1-(1/length(r1))
            r1 <- qnorm(r1/2+.5)
            r2 <- qnorm(r2/2+.5)
            dim(r1) <- dim(r2) <- dim(ss)
    		rownames(r1) <- rownames(r2) <- rownames(ss)
            pos <- match(names(tfmode), rownames(r1))
            r1 <- filterRowMatrix(r1, pos)
            r2 <- filterRowMatrix(r2, pos)
    		dnull <- sapply(1:per, function(i, ll, samp) {
    			ll[sample(length(ll), samp)] <- 0
    			return(ll)
    		}, ll=ll, samp=length(tgs))
    		dnull <- cbind(ll, dnull)
    		dnull[match(tgs, names(tfmode)), 1] <- 0
    		dnull <- t(t(dnull)/colSums(dnull))		
    		sum1 <- t(r2) %*% (dnull*tfmode)
            sum2 <- t(r1) %*% (dnull * (1-abs(tfmode)))
            es <- colMeans(abs(sum1)+sum2*(sum2>0))
    		x <- es[1]
    		es <- es[-1]
            iqr <- quantile(es, c(.5, 5/length(es)))
            pd <- ecdf(es)
            a <- list(x=knots(pd), y=pd(knots(pd)))
            filtro <- a$x<iqr[1] & a$x>=iqr[2] & a$y<1
            spl <- smooth.spline(a$x[filtro], -log(a$y[filtro]), spar=.75)
            p <- exp(-predict(spl, x)$y)
            pos <- which(x>iqr[1])
            if (x>iqr[1]) p <- pd(x)
            return(p)
    	}, ss=mobj$signature, nes=mobj$es$nes, regul=mobj$regulon, pos=pos, pb=pb, verbose=verbose)
    }
    if (verbose) message("\nProcess ended at ", date())
	names(res) <- names(mobj$regulon)[pos]
	mobj$es$synergy <- res
    class(mobj) <- "msviper"
	return(mobj)	
}


comregulationAnalysis <- function(tfs, ges, regulon, nullmodel=NULL, ones, resp, minsize=5, adaptive.size=FALSE, cores=1, processAll=FALSE, verbose=TRUE) {
    if (cores>1) {
        reg1 <- .viperParallelApply(1:ncol(tfs), function(i, tfs, regulon, res1) {
            x <- tfs[, i]
            pos <- which(names(regulon) %in% x)
            tgs <- table(unlist(lapply(regulon[pos], function(x) names(x$tfmode)), use.names=FALSE))
            tgs <- names(tgs)[tgs==length(x)]
            if (length(tgs)<2) return(list(tfmode=NULL, likelihood=NULL))
            tfmode <- t(t(sapply(regulon[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs))*sign(res1[pos]))
            likelihood <- apply(sapply(regulon[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs), 1, prod)
            list(tfmode=rowMeans(tfmode, na.rm=TRUE), likelihood=likelihood)
        }, tfs=tfs, regulon=regulon, res1=ones, mc.cores=cores)
    }
    else {
        reg1 <- apply(tfs, 2, function(x, regulon, res1) {
            pos <- which(names(regulon) %in% x)
            tgs <- table(unlist(lapply(regulon[pos], function(x) names(x$tfmode)), use.names=FALSE))
            tgs <- names(tgs)[tgs==length(x)]
            if (length(tgs)<2) return(list(tfmode=NULL, likelihood=NULL))
            tfmode <- t(t(sapply(regulon[pos], function(x, tgs) x$tfmode[match(tgs, names(x$tfmode))], tgs=tgs))*sign(res1[pos]))
            likelihood <- apply(sapply(regulon[pos], function(x, tgs) x$likelihood[match(tgs, names(x$tfmode))], tgs=tgs), 1, prod)
            list(tfmode=rowMeans(tfmode, na.rm=TRUE), likelihood=likelihood)
        }, regulon=regulon, res1=ones)
    }
    names(reg1) <- apply(tfs, 2, paste, collapse="--")
    reg1 <- reg1[sapply(reg1, function(x) length(x$tfmode))>0]
    if (adaptive.size) {
        reg1 <- reg1[sapply(reg1, function(x) {
            sum(x$likelihood/max(x$likelihood))
        })>=minsize]
    }
    else reg1 <- reg1[sapply(reg1, function(x) length(x$tfmode))>=minsize]
    if (length(reg1)==0) return(list())
    if (is.null(nullmodel)) nullf <- ppwea3NULLf(reg1)
    else nullf <- pwea3NULLf(pwea3NULLgroups(nullmodel, reg1, cores=cores, verbose=verbose), cores=cores, verbose=verbose)
    res1 <- groupPwea3(ges, reg1, nullf, minsize=1, cores=cores, verbose=verbose)
    filtro <- sapply(strsplit(names(res1$p.value), "--"), function(x, res1, res) {
        test <- res[sapply(strsplit(names(res), "--"), function(x1, x) all(x1 %in% x), x=x)]
        if (length(test)>0) return(all(test > res1[paste(x, collapse="--")]))
        return(FALSE)
    }, res1=res1$p.value, res=resp)
    if (processAll) filtro <- rep(TRUE, length(filtro)) # Force all TRUE to process all pairs regardeless of enrichment
    tmp <- lapply(res1, function(x, filtro) {
        ifelse(is.null(ncol(x)), return(x[filtro]), return(filterRowMatrix(x, filtro)))
    }, filtro=filtro)
    return(list(es=tmp, regul=reg1[names(reg1) %in% names(tmp$nes)]))
}

selectTFs <- function(res, n) {
	tfs <- unlist(lapply(strsplit(names(res$nes), "--"), function(x, tfs, n) {
		tfs <- unique(unlist(tfs[sapply(tfs, function(tfs, x) any(tfs %in% x), x=x)], use.names=FALSE))
		if (length(tfs)<n) return(NULL)
		return(combn(tfs, n))
	}, tfs=strsplit(names(res$nes), "--"), n=n), use.names=FALSE)
	tfs <- matrix(tfs, n, length(tfs)/n)
	tfs <- tfs[, !duplicated(apply(tfs, 2, paste, collapse="--"))]
	i <- 1
	while(i<ncol(tfs)) {
		filtro <- apply(tfs, 2, function(x, pat) all(x %in% pat), pat=tfs[, i])
		filtro[i] <- FALSE
		tfs <- filterColMatrix(tfs, !filtro)
		i <- i+1
	}
	return(tfs)
}

#' @method print msviper
#' @export
print.msviper <- function(x, ...) cat("Object of class msviper with ", length(x$es$nes), " regulators.\n", sep="")

#' List msviper results
#' 
#' This function generates a table of msviper results
#' 
#' @param object msviper object
#' @param mrs Either number of top MRs to report or vector containing the genes to display
#' @param ... Given for compatibility with the summary generic function
#' @return Data.frame with results
#' @method summary msviper
#' @export

summary.msviper <- function(object, mrs=10, ...) {
    if (length(mrs)==1) {
        mrs <- names(object$es$nes)[order(object$es$p.value)[1:round(mrs)]]
        mrs <- mrs[order(object$es$nes[match(mrs, names(object$es$nes))], decreasing=TRUE)]
    }
    pos <- match(mrs, names(object$es$nes))
    tmp <- data.frame(Regulon=mrs, Size=object$es$size[pos], NES=round(object$es$nes[pos], 2), p.value=signif(object$es$p.value[pos], 3), FDR=signif(p.adjust(object$es$p.value, "fdr")[pos], 3))
    if (!is.null(object$es$synergy)) {
        synp <- object$es$synergy[match(tmp$Regulon, names(object$es$synergy))]
        tmp <- cbind(tmp, Synergy=signif(synp, 3))
    }
    if (!is.null(object$ledge)) {
        tmp <- cbind(tmp, Ledge=sapply(object$ledge[match(tmp$Regulon, names(object$ledge))], function(x) {
            if (length(x)<6) return(paste(x, collapse=", "))
            return(paste(paste(x[1:4], collapse=", "), ", + ", length(x)-4, " genes", sep=""))
        }))
    }
    if (!is.null(object$shadow)) {
        tmps <- apply(object$shadow, 1, paste, collapse=" -> ")
        tmp <- list("msviper.results"=tmp, "Shadow.pairs"=tmps)
    }
    return(tmp)
}
