bgmyc.lik <-
function(params, data) {
			n<-data$n[[params[3]]]
			p <- c(rep(params[2], n), params[1])					
			lambda <- sum(data$list.s.nod[[params[3]]][1:n, ])/sum(data$list.i.mat[[params[3]]][1:n, ]^p[-(n + 1)] %*% data$internod)

			lambda <- c(rep(lambda, n), sum(data$list.s.nod[[params[3]]][n + 1, ])/data$list.i.mat[[params[3]]][n + 1, ]^p[n + 1] %*% data$internod) 
			b <- t(data$list.i.mat[[params[3]]]^p) %*% lambda			
			lik <- b * exp(-b * data$internod)					
			out<-sum(log(lik))
		if(is.nan(out)){
			return(-Inf)
			}
		else{
			return(sum(log(lik)))}
	}
