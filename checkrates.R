checkrates<-function(result){
  
  bgmyc.lambda<-function(params, data) {
    n<-data$n[[params[3]]]
    p <- c(rep(params[2], n), params[1])					
    lambda <- sum(data$list.s.nod[[params[3]]][1:n, ])/sum(data$list.i.mat[[params[3]]][1:n, ]^p[-(n + 1)] %*% data$internod)
    lambda <- c(lambda, sum(data$list.s.nod[[params[3]]][n + 1, ])/data$list.i.mat[[params[3]]][n + 1, ]^p[n + 1] %*% data$internod) 
    return(lambda[c(2,1)])
  }
  
  
  branchrates<-c()
  if(class(result)=="multibgmyc"){
    for(i in 1:length(result)){
      print(i)
      bgmyc.dataprep(result[[i]]$tree)->dat
      for(j in 1:length(result[[i]]$par[,1])){
        branchrates<-rbind(branchrates, c(result[[i]]$par[j,1:3], bgmyc.lambda(result[[i]]$par[j,1:3], dat)))
        } 
      }
    }
    if(class(result)=="singlebgmyc"){
      bgmyc.dataprep(result$tree)->dat
      for(j in 1:length(result$par[,1])){
        branchrates<-rbind(branchrates, c(result$par[j,1:3], bgmyc.lambda(result$par[j,1:3], dat)))
      }
    }
    
  ratemod<-((branchrates[,3]^branchrates[,1]))/branchrates[,3]
  branchrates<-cbind(branchrates, branchrates[,4]*ratemod) 
  
  colnames(branchrates)<-c("p.div", "p.coal", "threshold", "lambda.div", "lambda.coal", "lambda.div.mod")
  class(branchrates)<-"bgmycrates"
  return(branchrates)
  }