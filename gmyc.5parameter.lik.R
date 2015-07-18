gmyc.5parameter.lik <-
function(params, data) {
		p <- c(rep(params[2], data$n[[params[5]]]), params[1])					
		lambda <- c(rep(params[4], data$n[[params[5]]]), params[3]) 
		b <- t(data$list.i.mat[[params[5]]]^p) %*% lambda			
		lik <- b * exp(-b * data$internod)					
		return(sum(log(lik)))
	}
