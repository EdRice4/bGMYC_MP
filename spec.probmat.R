spec.probmat <-
function(res) {

	spec.list<-function(res, thresh) {
		tr <- res$tree
		spec <- tr$tip.label
		numtip <- length(tr$tip.label)
		

		max.mrca <- res$mrca[[thresh]] + numtip
		
		numspec <- length(max.mrca)
		
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
		
		res2 <- c()
		
		for (i in 1:length(max.mrca)) {
			tip.name <- tr$tip.label[nest.tip(max.mrca[i], tr)]
			res2 <- rbind(res2, cbind(i, tip.name))	
		}
		
		
		if (length(spec[-match(res2[,2], spec)]) != 0) {
			numspec <- numspec + 1
			for (s in spec[-match(res2[,2], spec)]) {
				res2 <- rbind(res2, cbind(numspec, s))
				numspec <- numspec + 1
			}
		}
		
		res2 <- data.frame(res2)
		colnames(res2) <- c("GMYC_spec", "sample_name")
		
		return (res2)
	}

assignlists<-c()

if(class(res)=="singlebgmyc"){
	numtip <- length(res$tree$tip.label)
	probmat<-matrix(data=( rep(0, times=(numtip*numtip))), ncol=numtip, nrow=numtip)
	rownames(probmat)<-res$tree$tip.label
	colnames(probmat)<-res$tree$tip.label

	for(j in 1:length(res$par[,3])){
		assignlists<-spec.list(res, res$par[j,3])
		levs<-levels(assignlists[,1])
		for(i in 1:length(levs)){
		probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]<-probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]+(1/length(res$par[,3]))
			}
		}
}

if(class(res)=="multibgmyc"){
	numtip <- length(res[[1]]$tree$tip.label)
	probmat<-matrix(data=( rep(0, times=(numtip*numtip))), ncol=numtip, nrow=numtip)
	rownames(probmat)<-res[[1]]$tree$tip.label
	colnames(probmat)<-res[[1]]$tree$tip.label
	ntrees<-length(res)
	totalsamp<-ntrees*length(res[[1]]$par[,1])

	for(q in 1:ntrees){
		for(j in 1:length(res[[q]]$par[,3])){
			assignlists<-spec.list(res[[q]], res[[q]]$par[j,3])
			levs<-levels(assignlists[,1])
			for(i in 1:length(levs)){
				probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]<-probmat[as.character(assignlists[which(assignlists[,1]==levs[i]),2]),as.character(assignlists[which(assignlists[,1]==levs[i]),2])]+1

				}
			}
	}
	probmat<-probmat/totalsamp
}

class(probmat)<-"bgmycprobmat"
return(probmat)

}
