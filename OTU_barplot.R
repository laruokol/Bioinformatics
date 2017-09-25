# OTU_barplot
# 
# This function plots the relative abundance of a taxon (e.g. OTU)
# as a barplot, where the samples are arranged according to a
# grouping factor.
#
# INPUT: physeq  a phyloseq-object
#        group   name of the grouping factor in the sample_data()
#        OTU.num row number in the otu_table() to be plotted
#
# (c) Lasse Ruokolainen -- November 2016
#####################################################################

OTU_barplot = function(physeq,group,OTU.num){
	require(RColorBrewer)
	Y = transform_sample_counts(physeq,function(x) x/sum(x))
	tmp = cbind(data.frame(sample_data(Y)),
				OTU=as.vector(otu_table(Y)[as.character(OTU.num)]))
	tmp = tmp[order(tmp[,group],tmp$OTU,decreasing=T),]
	cols = tmp[,group]
	if(nlevels(cols)==6){
		levels(cols) = c('lightgreen','forestgreen','lightgoldenrod',
						 'orange3','skyblue','royalblue')
	}else{
		levels(cols) = brewer.pal(nlevels(cols),'Paired')
	}	
	m = max(sqrt(tmp$OTU))
	
	par(xpd=NA)
	x=barplot(sqrt(tmp$OTU),col=as.character(cols),ylab='sqrt Relative abundance',
			names.arg=F,las=3,ylim=c(0,m),border=F)
	re = rle(as.numeric(tmp[,group]))
	s = cumsum(re$len)
	for(ii in 1:length(s)){
		lines(rep(mean(x[s[ii]:(s[ii]+1)]),2),c(-m/20,m),lwd=2,col='black')
	}
	for(ii in levels(tmp[,group])){
		w = tmp[,group]==ii
		text(mean(x[w]),-m/20,ii,font=3)
	}
}
