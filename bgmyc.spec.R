bgmyc.spec <-
function(res, filename=NULL, cmatrix=NULL) {
	temp.file.mcmc<-tempfile()
	
	if(class(res)=="singlebgmyc"){
		numsamp<-length(res$par[,1])
		tr <- res$tree
		spec <- tr$tip.label
		numtip <- length(tr$tip.label)
		output<-list()
		
		nest.tip <- function(nod, tr) {
			
			tip <- c()
			child <- tr$edge[tr$edge[,1] == nod, 2]
			
			for (ch in child) {
				if (ch <= numtip) {
					tip <- c(tip, ch)
				} else {
					tip <- c(tip, nest.tip(ch, tr))
				}
			}
			return (tip)
		}
		
		clusters<-list()

		for (i in 1:length(res$par[,3])){


			max.mrca <- res$mrca[[res$par[(i),3]]] + numtip
			numspec <- length(max.mrca)

			result <- c()
			delimit<-list()
			
			for(j in 1:numspec){
				tip.name <- tr$tip.label[nest.tip(max.mrca[j], tr)]
				result <- rbind(result, cbind(j, tip.name))
				cat(paste(unlist(sort(tip.name)), collapse="	"), "\n", file=temp.file.mcmc, append=TRUE)
				delimit[[j]]<-sort(tip.name)
			}

			if (length(spec[-match(result[,2], spec)]) != 0) {
				numspec <- numspec + 1	
				for (s in spec[-match(result[,2], spec)]) {
					result <- rbind(result, cbind(numspec, s))
					cat(s, "\n", file=temp.file.mcmc, append=TRUE)
					delimit[[numspec]]<-s
					numspec <- numspec + 1
					
				}
			}
			clusters[[i]]<-delimit
		}

		commandline<-c("sort ", temp.file.mcmc, " | uniq -c | sed s/'^ *'// | sed s/' $'//")
			
		read.table(pipe(paste(unlist(commandline), collapse = " ")), sep=" ")->output[["specprobs"]]
		closeAllConnections()

		output$specprobs<-output$specprobs[order(output$specprobs[,1], decreasing=TRUE),]
		output$specprobs[,1]<-output$specprobs[,1]/numsamp
		
		if(!is.null(filename)){
			write.table(output[["specprobs"]], filename, row.names=FALSE, col.names=FALSE)
			}
		
		
		matrices<-list()
		if(!is.null(cmatrix)){

			for(k in 1:length(clusters)){
				vec<-c()
				for(m in 1:length(clusters[[k]])){
					vec<-c(vec, clusters[[k]][[m]][1])
				}
				cols<-c(colnames(cmatrixname))

				newmat<-matrix(nrow=length(clusters[[k]]), ncol=ncol(cmatrixname), dimnames=(list(vec, cols) ))
				

				for(n in 1:length(clusters[[k]])){
					if(length(clusters[[k]][[n]])>1){			
						newmat[clusters[[k]][[n]][[1]],]<-colSums(cmatrixname[clusters[[k]][[n]],])	
					}else{newmat[clusters[[k]][[n]][[1]],]<-cmatrixname[clusters[[k]][[n]],]}
				}
				matrices[[k]]<-newmat
			}
			
		output[["matrices"]]<-matrices	
			
		}

	}

	if(class(res)=="multibgmyc"){
	
		ntrees<-length(res)	
		samp<-length(res[[1]]$par[,1])
		totalsamp<-length(res[[1]]$par[,1])*length(res)
		numtip <- length(res[[1]]$tr$tip.label)

		nest.tip <- function(nod, tr) {
			
			tip <- c()
			child <- tr$edge[tr$edge[,1] == nod, 2]
		
			for (ch in child) {
				if (ch <= numtip) {
					tip <- c(tip, ch)
				} else {
					tip <- c(tip, nest.tip(ch, tr))
				}
			}
			return (tip)
		}


		clusters<-list()
		for(h in 1:ntrees){
			tr <- res[[h]]$tree
			spec <- tr$tip.label
			output<-list()
		

			for (i in 1:samp){

				max.mrca <- res[[h]]$mrca[[res[[h]]$par[(i),3]]] + numtip
				numspec <- length(max.mrca)

				result <- c()
				delimit<-list()
			
				for(j in 1:numspec){
					tip.name <- tr$tip.label[nest.tip(max.mrca[j], tr)]
					result <- rbind(result, cbind(j, tip.name))
					cat(paste(unlist(sort(tip.name)), collapse="	"), "\n", file=temp.file.mcmc, append=TRUE)
					delimit[[j]]<-sort(tip.name)
				}

				if (length(spec[-match(result[,2], spec)]) != 0) {
					numspec <- numspec + 1	
					for (s in spec[-match(result[,2], spec)]) {
						result <- rbind(result, cbind(numspec, s))
						cat(s, "\n", file=temp.file.mcmc, append=TRUE)
						delimit[[numspec]]<-s
						numspec <- numspec + 1
					
					}
				}
				clusters[[i+(samp*(h-1))]]<-delimit
			}
		}

		matrices<-list()

		if(!is.null(cmatrix)){


			for(k in 1:length(clusters)){
				vec<-c()
				for(m in 1:length(clusters[[k]])){
					vec<-c(vec, clusters[[k]][[m]][1])
				}
				cols<-c(colnames(cmatrixname))

				newmat<-matrix(nrow=length(clusters[[k]]), ncol=ncol(cmatrixname), dimnames=(list(vec, cols) ))
				

				for(n in 1:length(clusters[[k]])){
					if(length(clusters[[k]][[n]])>1){			
						newmat[clusters[[k]][[n]][[1]],]<-colSums(cmatrixname[clusters[[k]][[n]],])	
					}else{newmat[clusters[[k]][[n]][[1]],]<-cmatrixname[clusters[[k]][[n]],]}
				}
				matrices[[k]]<-newmat
			}
			
		output[["matrices"]]<-matrices	
			
		}

		commandline<-c("sort ", temp.file.mcmc, " | uniq -c | sed s/'^ *'// | sed s/' $'//")
			
		read.table(pipe(paste(unlist(commandline), collapse = " ")), sep=" ")->output[["specprobs"]]
		closeAllConnections()

		output$specprobs<-output$specprobs[order(output$specprobs[,1], decreasing=TRUE),]
		output$specprobs[,1]<-output$specprobs[,1]/totalsamp
		
		if(!is.null(filename)){
			write.table(output[["specprobs"]], filename, row.names=FALSE, col.names=FALSE)
			}
			
	}
		unlink(temp.file.mcmc)
		return(output)
}
