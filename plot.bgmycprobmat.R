plot.bgmycprobmat <-
function(probmat, tree){
	
	ntip<-length(tree$tip.label)
	
	tree<-reorder.phylo(tree, order="cladewise")
	
	edge<-tree$edge[which(tree$edge[,2]<(ntip+1)),]
	tiplab<-tree$tip.label[edge[,2]]

	orderedmat<-probmat[tiplab,]
	orderedmat<-orderedmat[,tiplab]
	
	par(fig=c(0, 0.5, 0, 1))
	plot(tree, show.tip.label=FALSE, no.margin=TRUE)

	par(fig=c(0.47, 0.96, 0.035, 0.965), new=TRUE)
	image(x=c(1:ntip+1), y=c(1:ntip), z=orderedmat, axes=FALSE, breaks=c(0, 0.05, 0.49999999, 0.899999999, 0.9499999999, 1), col=heat.colors(5))
	
	par(fig=c(0.968, 0.99, 0.035, 0.965), new=TRUE)
	leg<-matrix(nrow=1, ncol=20, data=seq(from=0.025, to=1, by=0.05))
	image(y=seq(from=0, to=1, by=0.05), x=1, z=leg, col=heat.colors(5), breaks=c(0, 0.05, 0.5, 0.9, 0.95, 1), axes=FALSE)

	abline(h=c(0.95, 0.9, 0.5, 0.05))
	text(x=1, y=0.975, labels="p=\n0.95-1", cex=0.5, srt=90)
	text(x=1, y=0.925, labels="p=0.9-\n0.95", cex=0.5, srt=90)
	text(x=1, y=0.7, labels="p=0.5-0.9", cex=0.5, srt=90)
	text(x=1, y=0.275, labels="p=0.05-0.5", cex=0.5, srt=90)
	text(x=1, y=0.025, labels="p=0-\n0.05", cex=0.5, srt=90)
	}
