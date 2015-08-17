plot.multibgmyc <-
function(result, plot=TRUE){

	parmat<-c()
	for(i in 1:length(result)){
		parmat<-rbind(parmat, result[[i]]$par)
		}
		
	if(plot){
		par(mfrow=c(2,2))
		ylabels<-c("p.div","p.coal","threshold","logposterior")
		for(i in 1:4){
			plot(parmat[,i], xlab="generations", ylab=ylabels[i])
			}
		}
	
	else{return(parmat)}
}
