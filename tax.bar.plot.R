# tax.bar.plot
#
# This function plots relative abundances at a given taxonomy level.
#
# INPUT:    physeq    = phyloseq-class object
#           tax_level = desired taxonomy level (default = 'Phylum')
#           ncols     = how many most abundant taxa should be highlighted? 
#                       (default = 10)
#           ord       = optional grouping factor (samples will be arrangen 
#                       in a descending
#						order according to 'ord')
#           palette   = color palette to be used (uses predefined palettes 
#                       for brewer.pal)
#           names     = optional specification of sample names
#
# (c) Lasse Ruokolainen, April 2016
##########################################################################

  
tax.bar.plot = function(physeq, tax_level="Phylum", ncols=10, ord = NULL, palette='Paired',names=NULL,leg.pos.x=NULL,lmar=10) {
	require(phyloseq)
	require(RColorBrewer)
	tmp = physeq
	
	tax = as.factor(tax_table(physeq)[,tax_level])
	levels(tax) = ifelse(levels(tax)=='','Unknown',levels(tax))

	mat = matrix(0,nlevels(tax),nsamples(physeq))
	rownames(mat) = levels(tax)
	for(ii in levels(tax)){
		w = tax%in%ii
		if(sum(w)==1){
			mat[ii,] = otu_table(physeq)[w,]
		}else{
			mat[ii,] = colSums(otu_table(physeq)[w,])
		}		
	}
	y = apply(mat,2,function(x) x/sum(x))
	
	o = order(rowMeans(y),decreasing=T)

	if(is.null(ord)){
		oo = order(y[o[1],],decreasing=T)	
		gap = 0
	}else{
		oo = order(ord,y[o[1],],decreasing=T)
		tmp2 = rle(as.numeric(ord[oo]))
		gap = 0
		for(ii in 1:nlevels(ord)){
			gap = c(gap,rep(0,tmp2$len[ii]-1),1.5)
		}
		gap = gap[-length(gap)]
	}

	if(palette=='Paired'){
		palette = brewer.pal(12,palette)
		if(ncols<=12) tmp = ncols else tmp = 12
		palette = palette[c(1:4,7:10,5,6,11,12)][1:tmp]
	}
	
	getPalette = colorRampPalette(palette)
	cols = c(getPalette(ncols),rep('grey80',nrow(y)-5))
	y = rbind(y[o,oo][1:ncols,],colSums(y[o,oo][(ncols+1):nrow(y),]))
	
	if(!is.null(names)){
		if(names==F) names = rep('',ncol(y))	
	}		
	if(is.null(leg.pos.x)) leg.pos.x = max(x)*1.05

			
	op=par(mar=c(4,4,1,lmar),mgp=c(2,.8,0),xpd=NA,lwd=.5)
	x = barplot(y,col=cols,axes=F,las=2,cex.names=.5,
	            ylab='Proportional abundance',cex.lab=1.25,
	            space=gap,names.arg=names,border=NA)
	axis(2)
	legend(leg.pos.x,1,c(as.character(levels(tax)[o][1:ncols]),'Other'),
		   pch=15,col=cols,pt.cex=1.5, bty="n",title=tax_level,title.adj=0.05)		   

	if(!is.null(ord)){
		tmp = ord[oo]
		for(ii in levels(tmp)){
			w = tmp==ii
			text(mean(x[w]),-.05,ii,font=3)
		}
	}
	on.exit(par(op))
}
