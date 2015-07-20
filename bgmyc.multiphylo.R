bgmyc.multiphylo <- function(
                             multiphylo, mcmc, burnin, thinning, py1=0, py2=2,
                             pc1=0, pc2=2, t1=2, t2=51, scale=c(20, 10, 5.00),
                             start=c(1.0, 0.5, 50.0), sampler=bgmyc.gibbs,
                             likelihood=bgmyc.lik, prior=bgmyc.prior,
                             mpibool=FALSE)
{
	
    ntre <- length(multiphylo)

    cat("You are running a multi tree analysis on", ntre, "trees.\n")
    cat("These trees each contain", length(multiphylo$tip.label[[1]]), "tips.\n")
    cat("The Yule process rate change parameter has a uniform prior ranging from", py1, "to", py2, ".\n")
    cat("The coalescent process rate change parameter has a uniform prior ranging from", pc1, "to", pc2, ".\n")
    cat("The threshold parameter, which is equal to the number of species, has a uniform prior ranging from", t1, "to", t2, ". The upper bound of this prior should not be more than the number of tips in your trees.\n")
    cat("The MCMC will start with the Yule parameter set to", start[1], ".\n")
    cat("The MCMC will start with the coalescent parameter set to", start[2], ".\n")
    cat("The MCMC will start with the threshold parameter set to", start[3], ". If this number is greater than the number of tips in your tree, an error will result.\n")
    cat("Given your settings for mcmc, burnin and thinning, your analysis will result in", ((mcmc-burnin)/thinning)*ntre, "samples being retained.\n")

    # Test for MPI environment
    if (Sys.info()['sysname'] == "Linux") {
        nproc <- as.numeric(system("nproc", intern=T))
    }
    if (Sys.info()['sysname'] == "Darwin") {
        nproc <- as.numeric(system("sysctl -n hw.ncpu", intern=T))
    }

    if (mpibool && nproc >= 3) {
        # Spawn slave CPUs, preserving one for master
        Rmpi::mpi.spawn.Rslaves(nslaves=nproc-1)
        # Calculate how many trees to send to each (keeping portion for master)
        buffer <- ceiling(length(multiphylo) / nproc)
        # Partition data
        trees.split <- split(multiphylo, ceiling(seq_along(multiphylo) / buffer))
        # Scatter data among members in comm
        mpi.bcast.cmd(trees.ind <- mpi.scatter.Robj(trees.split))
        trees.ind <- mpi.scatter.Robj(multiphylosplit)

    outputlist <- list()
    
    for(i in 1:ntre){
        bgmyc.dataprep(multiphylo[[i]])->data
        NNodes<-data$tree$Nnode
        sampler(data, m=mcmc, burnin, thinning, py1=py1, py2, pc1, pc2, t1, t2, scale=scale, start=start, likelihood=likelihood, prior=prior)->outputlist[[i]]
        cat("tree number ", i, "is finished", "\n")
        }
    class(outputlist)<-"multibgmyc"
    return(outputlist)
}
