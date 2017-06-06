download.Yeast.Go.terms <-
function(url.GO.terms='http://downloads.yeastgenome.org/curation/literature/go_terms.tab'){
  #### input
  #### output 
  	
  GO.terms <- read.delim(url(url.GO.terms) 
                      ,header=FALSE,check.names=FALSE, stringsAsFactors=FALSE)
  colnames(GO.terms) <- c('GOID' ,'GO_Term', 'GO_Aspect', 'GO_Term_Definition')
  GO.terms
}
