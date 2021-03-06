#!/usr/bin/env Rscript

bgmyc.multiphylo.mpi <- function(
                                 multiphylo, mcmc, burnin, thinning, py1=0,
                                 py2=2, pc1=0, pc2=2, t1=2, t2=51,
                                 scale=c(20, 10, 5), start=c(1, 0.5, 50),
                                 sampler=bgmyc.gibbs.mpi, likelihood=bgmyc.lik,
                                 prior=bgmyc.prior
                                 ) {

    # Test for MPI environment and determine number of CPUs to utilize
    #if (Sys.info()['sysname'] == "Linux") {
        #nproc <- as.numeric(system("nproc", intern=T))
    #}
    #if (Sys.info()['sysname'] == "Darwin") {
        #nproc <- as.numeric(system("sysctl -n hw.ncpu", intern=T))
    #}

    # Halt execution of script if insufficient amount of CPUs
    if (nproc <= 2) {
        stop(
             "This system has an insufficient number of CPUs.
             If running on linux, check no. of CPUs in terminal with 'nproc'.
             If running on mac, check no. of CPUs in terminal with
             'sysctl -n hw.cpu'.\n"
             )
    }

    # Calculate how many trees to send to each slave
    buffer <- ceiling(length(multiphylo) / (nproc - 1))
    # Partition data
    ntre <- length(multiphylo)
    trees.split <- split(multiphylo, ceiling(seq_along(multiphylo) / buffer))

    # Print informative output for user
    cat("You are running a multi tree analysis on", ntre, "trees.\n")
    cat("These trees each contain", length(multiphylo$tip.label[[1]]), "tips.\n")
    cat("The Yule process rate change parameter has a uniform prior ranging from", py1, "to", py2, ".\n")
    cat("The coalescent process rate change parameter has a uniform prior ranging from", pc1, "to", pc2, ".\n")
    cat("The threshold parameter, which is equal to the number of species, has a uniform prior ranging from", t1, "to", t2, ". The upper bound of this prior should not be more than the number of tips in your trees.\n")
    cat("The MCMC will start with the Yule parameter set to", start[1], ".\n")
    cat("The MCMC will start with the coalescent parameter set to", start[2], ".\n")
    cat("The MCMC will start with the threshold parameter set to", start[3], ". If this number is greater than the number of tips in your tree, an error will result.\n")
    cat("Given your settings for MCMC, burnin and thinning, your analysis will result in", ((mcmc-burnin)/thinning)*ntre, "samples being retained.\n")
    cat("Given your settings for MPI, your analysis will result in ")
    for(i in 1:length(trees.split)) {
        cat(length(trees.split[i]), "samples being sent to slave", i, ", ")
    }

    # Spawn slave CPUs, preserving one for master
    mpi.spawn.Rslaves(nslaves=nproc-1)

    # Optimize function for MPI environment
    bgmyc.mpi <- function(trees, ...) {

        # Initialize empty list for output
        outputlist <- list()

        for (i in 1:length(trees)) {
            data <- bgmyc.dataprep(trees[[i]])
            NNodes <- data$tree$Nnode
            outputlist[[i]] <- sampler(
                                       data, m=mcmc, burnin, thinning, py1,
                                       py2, pc1, pc2, t1, t2, scale, start,
                                       likelihood, prior
                                       )
        }

        class(outputlist) <- "multibgmyc"
        return(outputlist)

    }

    # Prepare environment; this is ugly
    mpi.bcast.Robj2slave(bgmyc.gibbs.mpi)
    mpi.bcast.Robj2slave(bgmyc.lik)
    mpi.bcast.Robj2slave(bgmyc.prior)
    mpi.bcast.Robj2slave(bgmyc.dataprep)
    mpi.bcast.Robj2slave(is.ultrametric)
    mpi.bcast.Robj2slave(is.binary.tree)
    mpi.bcast.Robj2slave(branching.times)
    
    # Run that function, boi
    output <- mpi.apply(
                        trees.split, bgmyc.mpi, mcmc, burnin, thinning, py1,
                        py2, pc1, pc2, t1, t2, scale, start, likelihood,
                        prior
                        )
    
    # Exit MPI
    mpi.exit()

    # Collapse output by first level only
    output <- unlist(output, recursive=FALSE)

    # Return output
    return(output)
}
