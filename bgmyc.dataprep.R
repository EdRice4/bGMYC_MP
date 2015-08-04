bgmyc.dataprep <- function(tr)
{

	if (!is.ultrametric(tr)) {
		stop("Your input tree is not ultrametric. This method requires that trees be ultrametric.")
	}
	if (!is.binary.tree(tr)) {
		stop("Your input tree is not fully bifurcating, please resolve with zero branch lengths")
	}
	if (0 %in% tr$edge.length[which(tr$edge[,2]<=length(tr$tip.label))]) {
		stop("Your tree contains tip branches with zero length. This will wreak havoc with the GMYC model.")
	}


	local.env <- environment()
	read.data <- function(z = 1) {
		bt <- -branching.times(tr)
		bt[bt > -1e-06] <- -1e-06
		names(bt) <- NULL
		assign("bt", bt, envir = local.env)
		assign("sb", sort(bt), envir = local.env)
		assign("numnod", length(bt), envir = local.env)
		assign("numtip", length(tr$tip.label), envir = local.env)
		assign("numall", length(bt) + length(tr$tip.label), envir = local.env)
		assign("nthresh", numnod, envir = local.env)
		
		internod <- sb[2:numnod] - sb[1:numnod - 1]
		internod[numnod] <- 0 - sb[numnod]
		assign("internod", internod, envir = local.env)

		assign("nesting", sapply((numtip + 1):numall, nesting.nodes), 
			envir = local.env)
		assign("nested", sapply((numtip + 1):numall, nest.nodes), 
			envir = local.env)

		ancs <- cbind(tr$edge[pmatch((1:numnod + numtip), tr$edge[, 
			2]), 1], (1:numnod + numtip))
		bt.ancs <- cbind(bt[ancs[, 1] - numtip], bt[ancs[, 2] - 
			numtip])
		assign("bt.ancs", bt.ancs, envir = local.env)
		
	}

					  nest.nodes <- function(x, p = 0) {
						  numtip <- length(tr$tip.label)
						  nods <- array(NA, 0)
						  desc <- as.integer(tr$edge[, 2][tr$edge[, 1] == x])
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
						  numtip <- length(tr$tip.label)
						  nod <- array(NA, 0)
						  if (x >= numtip + 2) {
							  anc <- as.integer(tr$edge[, 1][tr$edge[, 2] == x])
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
			threshy <- sb[j]									
			tmp <- (bt.ancs[, 1] < threshy) & (bt.ancs[, 2] >= threshy)		
			nod.type <- tmp + (bt >= threshy)								
			mrca.nodes[[j]] <- which(nod.type == 2)				
			cat(nod.type[1])
			if (nod.type[2] == 1) {
				nod.type[1] <- 2
			}
			nod.types[[j]]<-nod.type		
			
	   		n[[j]]<-length(mrca.nodes[[j]])		
	   		
			list.i.mat[[j]] <- matrix(0, ncol = numnod, nrow = (n[[j]] + 1))			
			list.s.nod[[j]] <- matrix(0, ncol = numnod, nrow = (n[[j]] + 1))			

			nod[[j]]<-nod.types[[j]][order(bt)]	

			for (i in (1:n[[j]])) {										
				list.s.nod[[j]][i, mrca.nodes[[j]][i]] <- 2										
 				if (!is.null(nested[[mrca.nodes[[j]][i]]])) {				
					list.s.nod[[j]][i, nested[[mrca.nodes[[j]][i]]] - numtip] <- 1		
				}
				list.s.nod[[j]][i, ] <- list.s.nod[[j]][i, order(bt)]		
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
		prepdata[["tree"]]<-tr

	
	return(prepdata)
}
