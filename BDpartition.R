# BDpartition
#
# a function to partition beta diversity in a  
# community data. Parameter n gives the number 
# of taxa and samples plotted (default = c(50,30)).
#
# The data should be transformed before analysis 
# (default = 'hellinger'), see help for 'decostand' 
# in 'vegan'. For methodological details, see:
# Legendre & De Caceres 2013, Ecology Letters
#                                            
# (c) Lasse Ruokolainen 2017
####################################################


BDpartition = function(data,n=NULL,transform='hellinger'){ 
	library(phyloseq)
	library(vegan)
	
	if(class(data)=='phyloseq'){
		X = t(otu_table(physeq))
	}else{
		if(class(data)=='matrix' | class(data)=='data.frame'){
			X = data
		}else{
			error('data must be either a matrix, data.frame, 
			or a phyloseq object')
		}
	}
	
	if(is.null(n)){
		n = numeric(2)
		if(ncol(X)<50) n[1] = ncol(X) 
		if(nrow(X)<30) n[2] = nrow(X) 
	}
	
	Y = decostand(X,transform)
	S = apply(Y,2,function(x) (x-mean(x))^2)
	SStotal = sum(S)
	
	BDtotal = SStotal / (nrow(Y)-1)
	
	layout(matrix(c(1,2),1,2))
	
	# Species contribution:
	SCBDj = 100*colSums(S) / SStotal
	o = order(SCBDj,decreasing = F)
	op=par(mar=c(6,6,1,1),mgp=c(3,.6,0),cex.axis=.7,cex.lab=.8)
	barplot(SCBDj[o][(length(SCBDj)-n[1]+1):length(SCBDj)],
		horiz = T,las=2,col='wheat',
	        xlab=expression(paste('Species contribution to ',
	        beta,'-diversity (%)')))
	lines(rep(mean(SCBDj),2),c(0,n[1]+2*n[1]/10),col='firebrick',lwd=3)
	par(op)
	
	# Sample contribution:
	LCBDi = 100*rowSums(S) / SStotal
	o = order(LCBDi,decreasing = F)
	op=par(mar=c(6,6,1,1),mgp=c(3,.6,0),cex.axis=.7,cex.lab=.8)
	barplot(LCBDi[o][(length(LCBDi)-n[2]+1):length(LCBDi)],
	horiz = T,las=2,col='wheat',
	        xlab=expression(paste('Sample contribution to ',
	        beta,'-diversity (%)')))
	lines(rep(mean(LCBDi),2),c(0,n[2]+2*n[2]/10),col='firebrick',lwd=3)
	par(op)
	
	BDpartition = list(BDtotal=BDtotal, SCBDj=SCBDj,LCBDi=LCBDi)
	return(BDpartition)
}