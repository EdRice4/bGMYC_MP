bgmyc.prior <-
function(params, py1, py2, pc1, pc2, t1, t2){
	loglik<-dunif(params[1], min=py1, max=py2, log=TRUE)+
			dunif(params[2], min=pc1, max=pc2, log=TRUE)+
			dunif(params[3], min=t1, max=t2, log=TRUE)
	return(loglik)
	}
