# tax_glom_fast
#
# The function takes in a phyloseq-object
# and returns a matrix with with pooled abundances at 
# a desired taxonomic level (default = Phylum). The 
# default is to sum over OTUs, but other functions
# can also be used.
#
# (c) lasse Ruokolainen 2016
# last modified: 28.11.2017
#####################################################

tax_glom_fast = function(physeq,tax_level='Phylum',FUN='sum'){
	tax = as.factor(tax_table(physeq)[,tax_level])
	
	levels(tax) = ifelse(levels(tax)=='','Unknown',levels(tax))

	g = list(tax)
	mat = aggregate(otu_table(physeq),by=g,FUN)
	rownames(mat) = mat[,1]; mat = mat[,-1]
	return(mat)	
}
