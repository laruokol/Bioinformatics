# diversity.rarefaction.curves
#
# This function perfomrs a rarefaction of OTU data, within
# a phyloseq-object, and calculates (by default) species
# richness (q = 0), Shannon entropy (q = 1), and inverse
# Simpson index (q = 2). The calculation is based on Hill's
# numbners (q), which can be used to give an arbitrary weight
# to species relative abundaces in the calculation.
#
# INPUT:    phyloseq    a phyloseq-object with an OTU-table
#           qvals       Hill-parameters used
#           resolution  how many rarefaction steps to use?
#                       The minimum is 1/10th of min(sample 
#                       size) and max is min(sample size).
#           plot        should the results be plotted?
#                       (default = TRUE).
#           values      should the values be returned?
#                       (default = TRUE).
#
# (c) Lasse Ruokolainen -- October 2017
##############################################################

diversity.rarefaction.curves = function(phyloseq,qvals = c(0,1,2),resolution = NULL,plot=TRUE,values=TRUE){		
	# load required packages:
	library(phyloseq)
	library(vegan)
	library(tidyr)
	library(dplyr)
	library(ggplot2)
	library(foreach)
	library(doMC)
	registerDoMC(detectCores(logical=F))
	
	Y = phyloseq; n = nsamples(Y)
	if(is.null(resolution)) res = 10 else res = resolution

	tmp = signif(min(sample_sums(Y)),1)
	tmp = floor(min(sample_sums(Y))/1000)*1000
	rvals = seq(tmp/10,tmp,len = res)
	
	divs = foreach(ii = 1:res) %dopar% {
		# rarefy to even depth:
		X = rrarefy(t(otu_table(Y)),rvals[ii]) 
		
		tmp = matrix(0,length(qvals),n)
		for(jj in 1:length(qvals)){
			# calculate diversity of order q:
			tmp[jj,] = true.diversity(X,q = qvals[jj]) 
		}
		return(tmp)
	}
	
	D = matrix(0,n*res,length(qvals))
	colnames(D) = qvals
	for(ii in 1:ncol(D)){
		D[,ii] = matrix(sapply(divs,function(x) x[ii,]),n*res,byr=T)	
	}
	X = data.frame(rarv=rep(rvals,each=n),D,check.names=F)
	X = gather(X,key=rarv,value='Diversity'); names(X)[2] = 'q_value'
	
	if(values==TRUE) return(X)

	cc = c('salmon','olivedrab','skyblue')

	if(plot==TRUE){
	ggplot(X,aes(rarv,Diversity,color=q_value)) +
		geom_smooth() + scale_color_manual(values=cc,name='Hill parameter') +
		theme_bw() + theme(legend.position = c(0.1,.9),legend.background=element_blank(),
		legend.key=element_blank(),axis.title=element_text(size=14)) + 
		guides(color=guide_legend(override.aes=list(fill=NA))) +
		xlab('Number of reads sampled')
	}
}