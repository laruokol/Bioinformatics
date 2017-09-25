# tax_glom_fast
#
# The function takes in a phyloseq-object
# and returns a matrix with with pooled abundances at 
# a desired taxonomic level (default = Phylum).
#
# (c) lasse Ruokolainen 2016
#####################################################

tax_glom_fast = function(physeq,tax_level='Phylum'){
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
	colnames(mat) = sample_names(physeq)
	return(mat)	
}
