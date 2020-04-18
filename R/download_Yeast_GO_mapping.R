download_Yeast_GO_mapping <-
function(yeast.GO.url='http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'){
  ##### download Yeast Gene Ontology mapping  http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz
  GO.assocs.all  <- read.table(textConnection(readLines(gzcon(url(yeast.GO.url))))
                         ,header=FALSE,check.names=FALSE, stringsAsFactors=FALSE, sep="\t", quote=NULL, comment.char ='!'
                         ,col.names=c("Database","SGDID","DB_Object_Symbol","NOT","GOID","DB:Ref","Evidence"
                                     ,"With:From","GO_Aspect","DB_Object_Name","Synonym","Type","Taxon","Date","Asigned By","Notes 1","Notes 2"))
  ### use only most important columns "SGDID", "DB_Object_Symbol", "GOID"
  GO.assocs <- GO.assocs.all[,match(c("SGDID", "GOID"),colnames(GO.assocs.all))]
  ### one can filter for a certain evidence like IDA : Inferred from Direct Assay 
  ### in our case we will accept any of the Evidences and keep unique records of the association
  GO.assocs  <- unique(GO.assocs) #### for some reasons we have same rows duplications (for example the multiple Evidence levels: IBA, IC, IDA etc)
}
