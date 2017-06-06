convert.Yeast.SGDId.2.systematic <-
function(SGD.Ids=c("S000007287", "S000007287", "S000007287")){
  #### optional add the conversion from SGD Id to systematic Yeast gene names Y.... using synergizer (http://llama.mshri.on.ca/synergizer/translate/)
  ####  {"method":"translate","params":[{"authority":"ensembl","species":"Homo sapiens","domain":"hgnc_symbol","range":"entrezgene","ids":["snph","chac1","actn3","maybe_a_typo","pja1","prkdc","RAD21L1","Rorc","kcnk16"]}],"id":0}
  ####  input ....
  ####  output ....
  #### performs translation of iD from to using authority ... implemented by Synergizer (url...)
  
  
  ### this works OK:# curl -H "Content-Type: application/json" -d '{"method":"version", "params":[], "id":0}'  http://llama.mshri.on.ca/cgi/synergizer/serv
  url <- 'http://llama.mshri.on.ca/cgi/synergizer/serv'
  body <- '{"method":"version", "params":[], "id":0}';
  response<-POST(url, body = body , encode = "json", content_type_json())


  body <- paste('{"method":"translate","params":[{"authority":"sgd","species":"Saccharomyces cerevisiae","domain":"sgdid","range":"systematic","ids":["'
  		       ,paste(SGD.Ids,collapse='","')
  		       ,'"]}],"id":0}',sep='')
  response<-POST(url, body = body , encode = "json", content_type_json())
  if( !is.null(fromJSON(content(response,as="text"))$error) ){
  	stop(paste('ERROR when calling translation by Synergizer :',fromJSON(content(response,as="text"))$error))
  } else {
  	converted <- fromJSON(content(response,as="text"))$result
  	#colnames(converted) <- c("SGD.Id","Systematic")
  	return(converted)
  }		

}
