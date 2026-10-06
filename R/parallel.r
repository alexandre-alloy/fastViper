# Cross-platform replacement for the package's previous mclapply calls.
# PSOCK workers receive only the free variables used by FUN, rather than the
# entire calling frame. The active option prevents nested calls from starting
# additional worker pools.
.viperParallelNamespace <- function(env) {
    while (is.environment(env) && !identical(env, emptyenv())) {
        if (isNamespace(env)) return(getNamespaceName(env))
        env <- parent.env(env)
    }
    NULL
}

.viperParallelPrepare <- function(FUN) {
    oldenv <- environment(FUN)
    pkg <- .viperParallelNamespace(environment(.viperParallelPrepare))
    source.mode <- is.null(pkg)
    ns <- if (source.mode) NULL else asNamespace(pkg)
    worker.env <- new.env(parent=if (source.mode) baseenv() else ns)
    captures <- list()
    packages <- if (source.mode) character() else pkg
    queue <- list(list(fun=FUN, env=oldenv))
    while (length(queue)) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        globals <- codetools::findGlobals(current$fun, merge=FALSE)
        wanted <- unique(c(globals$variables, globals$functions))
        for (name in wanted) {
            if (exists(name, envir=worker.env, inherits=FALSE)) next
            env <- current$env
            while (is.environment(env) && !identical(env, emptyenv())) {
                if (!source.mode && identical(env, ns)) break
                if (source.mode && identical(env, baseenv())) break
                if (exists(name, envir=env, inherits=FALSE)) {
                    value <- get(name, envir=env, inherits=FALSE)
                    if (is.function(value) && !is.primitive(value)) {
                        fenv <- environment(value)
                        fpackage <- .viperParallelNamespace(fenv)
                        if (is.null(fpackage) && source.mode) {
                            cloned <- value
                            environment(cloned) <- worker.env
                            assign(name, cloned, envir=worker.env)
                            captures[[name]] <- cloned
                            queue[[length(queue) + 1L]] <- list(fun=value, env=fenv)
                        } else {
                            assign(name, value, envir=worker.env)
                            captures[[name]] <- value
                            if (!is.null(fpackage)) packages <- unique(c(packages, fpackage))
                        }
                    } else {
                        assign(name, value, envir=worker.env)
                        captures[[name]] <- value
                    }
                    break
                }
                env <- parent.env(env)
            }
        }
    }
    environment(FUN) <- worker.env
    list(fun=FUN, packages=packages, captures=captures)
}

.viperParallelApply <- function(tasks, FUN, ..., mc.cores=1L,
                                mc.preschedule=TRUE, mc.set.seed=TRUE) {
    n <- length(tasks)
    cores <- suppressWarnings(as.integer(mc.cores[1L]))
    if (is.na(cores) || cores < 1L) cores <- 1L
    cores <- min(cores, n)
    requested.cores <- cores
    if (!n) return(list())
    if (cores <= 1L || isTRUE(getOption("viper.parallel.active", FALSE)))
        return(lapply(tasks, FUN, ...))

    if (.Platform$OS.type != "windows") {
        old.option <- getOption("viper.parallel.active", FALSE)
        options(viper.parallel.active=TRUE)
        on.exit(options(viper.parallel.active=old.option), add=TRUE)
        return(parallel::mclapply(tasks, FUN, ..., mc.cores=cores,
                                  mc.preschedule=mc.preschedule,
                                  mc.set.seed=mc.set.seed))
    }

    prepared <- .viperParallelPrepare(FUN)
    args <- list(...)
    input.bytes <- as.numeric(utils::object.size(args)) +
        as.numeric(utils::object.size(prepared$captures))
    task.bytes <- as.numeric(utils::object.size(tasks))
    max.bytes <- getOption("viper.max.parallel.memory", 64 * 1024^3)
    if (length(max.bytes) != 1L || is.na(max.bytes) || max.bytes <= 0)
        max.bytes <- Inf
    # Approximate each worker's resident input copy. This is deliberately
    # conservative; users can raise the option on large-memory systems.
    fit <- function(k) input.bytes + task.bytes / k
    while (cores > 1L && cores * fit(cores) > max.bytes) cores <- cores - 1L
    if (cores < requested.cores)
        warning("Reducing PSOCK workers to ", cores, " to stay within the estimated viper.max.parallel.memory budget.", call.=FALSE)
    if (cores <= 1L) return(lapply(tasks, FUN, ...))

    old.option <- getOption("viper.parallel.active", FALSE)
    options(viper.parallel.active=TRUE)
    on.exit(options(viper.parallel.active=old.option), add=TRUE)
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add=TRUE)
    init.worker <- function(libpaths, pkgs) {
        .libPaths(libpaths)
        for (pkg in pkgs) loadNamespace(pkg)
        options(viper.parallel.active=TRUE)
        NULL
    }
    environment(init.worker) <- baseenv()
    parallel::clusterCall(cl, init.worker, .libPaths(), prepared$packages)
    if (mc.set.seed) parallel::clusterSetRNGStream(cl)
    task.worker <- function(x, worker.fun, args) {
        old <- getOption("viper.parallel.active", FALSE)
        options(viper.parallel.active=TRUE)
        on.exit(options(viper.parallel.active=old), add=TRUE)
        do.call(worker.fun, c(list(x), args))
    }
    environment(task.worker) <- baseenv()
    if (mc.preschedule)
        parallel::parLapply(cl, tasks, task.worker, worker.fun=prepared$fun, args=args)
    else
        parallel::parLapplyLB(cl, tasks, task.worker, worker.fun=prepared$fun, args=args)
}
