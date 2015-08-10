#!/usr/bin/env Rscript

# {{{  bgmyc.dataprep
bgmyc.dataprep <- function(tree)
{

	# {{{ Error
	# If the tree is not ultramteric, halt execution and display appropriate
	# error message.
	if (!is.ultrametric(tree)) {
		stop(
		     "Your input tree is not ultrametric. This method requires that
		     trees be ultrametric."
		     )
	}
	# If the tree is not bifurcating, halt execution and display appropriate
	# error message.
	if (!is.binary.tree(tree)) {
		stop(
		     "Your input tree is not fully bifurcating, please resolve with
		     zero branch lengths"
		     )
	}
	# If zero branch lengths in tree, halt execution and display appropriate
	# error message.
	if (0 %in% tree$edge.length[which(tree$edge[,2]<=length(tree$tip.label))]) {
		stop(
		     "Your tree contains tip branches with zero length. This will
		     wreak havoc with the GMYC model."
		     )
	}
	# }}}


	local.env <- environment()
	read.data <- function() {
		branch.times <- -branching.times(tree)
		branch.times[branch.times > -1e-06] <- -1e-06
		names(branch.times) <- NULL
		assign("branch.times", branch.times, envir = local.env)
		assign("sorted.branch.times", sort(branch.times), envir = local.env)
		assign("numnod", length(branch.times), envir = local.env)
		assign("numtip", length(tree$tip.label), envir = local.env)
		assign("numall", length(branch.times) + length(tree$tip.label), envir = local.env)
		assign("nthresh", numnod, envir = local.env)
		
		internod <- sorted.branch.times[2:numnod] - sorted.branch.times[1:numnod - 1]
		internod[numnod] <- 0 - sorted.branch.times[numnod]
		assign("internod", internod, envir = local.env)

		assign("nesting", sapply((numtip + 1):numall, nesting.nodes), 
			envir = local.env)
		assign("nested", sapply((numtip + 1):numall, nest.nodes), 
			envir = local.env)

		ancs <- cbind(
		              # Get corresponding beginning node for each end node
		              tree$edge[
		                      # Return position of each end node
		                      pmatch(
		                             # Labelling nodes
		                             (1:numnod + numtip),
		                             # End nodes
		                             tree$edge[, 2]
		                             )
		                      , 1]
		              , (1:numnod + numtip))
		# Issue is in pmatch
		## Issue is in the manner in which the nodes are labelled; the first
		## node of the tree, that posessing the lowest value,
		## the one furthest left, will never? be an end node and will therefore
		## never be present in the 2nd column of tree$edge.
		# a <- 1:numnod
		# b <- a + numtip
		branch.times.ancs <- cbind(branch.times[ancs[, 1] - numtip], branch.times[ancs[, 2] - 
			numtip])
		assign("branch.times.ancs", branch.times.ancs, envir = local.env)
		
	}

					  nest.nodes <- function(x, p = 0) {
						  numtip <- length(tree$tip.label)
						  nods <- array(NA, 0)
						  desc <- as.integer(tree$edge[, 2][tree$edge[, 1] == x])
						  if (desc[1] > numtip) {
							  nods <- c(nods, desc[1], nest.nodes(desc[1]))
						  }
						  if (desc[2] > numtip) {
							  nods <- c(nods, desc[2], nest.nodes(desc[2]))
						  }
						  if (length(nods) > 0) {
							  return(nods)
						  }
						  else {
							  return(NULL)
						  }
					  }

					  nesting.nodes <- function(x, p = 0) {
						  numtip <- length(tree$tip.label)
						  nod <- array(NA, 0)
						  if (x >= numtip + 2) {
							  anc <- as.integer(tree$edge[, 1][tree$edge[, 2] == x])
						  }
						  else {
							  anc <- 1
						  }
						  if (anc >= numtip + 2) {
							  nod <- c(nod, anc, nesting.nodes(anc))
						  }
						  if (anc == numtip + 1) {
							  nod <- c(nod, anc)
						  }
						  if (length(nod) > 0) {
							  return(nod)
						  }
						  else {
							  return(NULL)
						  }
					  }


	create.mat <- function(nthresh=numnod)
	{
	
		mrca.nodes<-list()
		nod.types<-list()
		n<-list()
		list.i.mat<-list()
		list.s.nod<-list()
		nod<-list()
		
		for (j in (2:nthresh)) {										
			# Threshy is the distinction?
			threshy <- sorted.branch.times[j]									
			# Tmp does not care about NA
			tmp <- (branch.times.ancs[, 1] < threshy) & (branch.times.ancs[, 2] >= threshy)		
			nod.type <- tmp + (branch.times >= threshy)								
			mrca.nodes[[j]] <- which(nod.type == 2)				
            if (nod.type[1] == 1) {
			#if (nod.type[2] == 1) {
				nod.type[1] <- 2
			}
			nod.types[[j]]<-nod.type		
			
	   		n[[j]]<-length(mrca.nodes[[j]])		
	   		
			list.i.mat[[j]] <- matrix(0, ncol = numnod, nrow = (n[[j]] + 1))			
			list.s.nod[[j]] <- matrix(0, ncol = numnod, nrow = (n[[j]] + 1))			

			nod[[j]]<-nod.types[[j]][order(branch.times)]	

			for (i in (1:n[[j]])) {										
				list.s.nod[[j]][i, mrca.nodes[[j]][i]] <- 2										
 				if (!is.null(nested[[mrca.nodes[[j]][i]]])) {				
					list.s.nod[[j]][i, nested[[mrca.nodes[[j]][i]]] - numtip] <- 1		
				}
				list.s.nod[[j]][i, ] <- list.s.nod[[j]][i, order(branch.times)]		
				list.i.mat[[j]][i, ][list.s.nod[[j]][i, ] == 2] <- 2		
				list.i.mat[[j]][i, ][list.s.nod[[j]][i, ] == 1] <- 1		
				list.i.mat[[j]][i, ] <- cumsum(list.i.mat[[j]][i, ])		
			}
			list.s.nod[[j]][list.s.nod[[j]] == 2] <- 1								

			list.i.mat[[j]] <- list.i.mat[[j]] * (list.i.mat[[j]] - 1)					
			list.s.nod[[j]][n[[j]] + 1, ] <- nod[[j]] == 0
			list.i.mat[[j]][n[[j]] + 1, nod[[j]] == 0] <- 1
			list.i.mat[[j]][n[[j]] + 1, nod[[j]] == 2] <- -1
			list.i.mat[[j]][n[[j]] + 1, ] <- cumsum(list.i.mat[[j]][n[[j]] + 1, ]) + 1
			
			

			
		}
		
		assign("mrca.nodes", mrca.nodes, envir = local.env)	
		assign("nod.types", nod.types, envir = local.env)
		assign("n", n, envir=local.env)
		assign("list.s.nod", list.s.nod, envir = local.env)
		assign("list.i.mat", list.i.mat, envir = local.env)
	}
	read.data()
	create.mat()
	
		prepdata<-list()
		prepdata[["mrca.nodes"]]<-mrca.nodes
		prepdata[["nod.types"]]<-nod.types
		prepdata[["n"]]<-n
		prepdata[["list.s.nod"]]<-list.s.nod
		prepdata[["list.i.mat"]]<-list.i.mat
		prepdata[["internod"]]<-internod
		prepdata[["tree"]]<-tree

	
	return(prepdata)
}
# }}}
