attribute.mut.inf <-
function(entity.attribute, show.progress=FALSE, alternative.calc=FALSE){
    ne <- length(unique(entity.attribute[,1])) ### number of unique entities    
    na <- length(unique(entity.attribute[,2])) ### number of unique attributes    

## it's faster to calculate entropy directly from definition since we have the frequency of each entities (1 for entity present 0 for absent)
#system.time({
  s<- tapply(entity.attribute[,1],entity.attribute[,2],length)/ne  ### this gives the probaility (frequency) to find any entity for each attribute
      ha <- -s*log(s) - (1-s)*log(1-s) ### entropy for each independent attribute (assumed to be on columns)
    #})



### if we have already the mutual information
    mi <- matrix(NA, ncol=length(ha), nrow=length(ha)); colnames(mi)<-rownames(mi)<-names(ha)

##  entropy between attributes and mutual information (can be parallelized)
total.iter <- (length(ha)-1) * length(ha)/2
 if(alternative.calc){
      ### this version is trying to avoid use of infotheo library (no need for calling entropy() function)
  for(j1 in 1:(length(ha)-1)){
    if( show.progress & ((j1 %% 10) == 0)) print(sprintf('completed:%.2f %%',j1*(j1-1)/2/total.iter*100))
    a1 <- names(ha)[j1]
    #system.time(
    for(j2 in (j1+1):length(ha)){
      a2 <- names(ha)[j2]
  	  ### get the mutual info between the 2 attributes
  	  eat <- table(entity.attribute[entity.attribute[,2]==a1 | entity.attribute[,2]==a2,])
  	  p <- table((eat[,1]*2+eat[,2]))/ne
  	  p <- c(1-sum(p),p)
  	  ### most of  the p are the same (most of the mi are the same for the same j1 !!!)	  	  
  	  mi[j1,j2] <- ha[j1]+ha[j2]- sum(-p* log(p))
  	}
    #)	
  }
  if( show.progress ) print(sprintf('completed:%.2f %%',100))
} else {	
  ent.att <- table(entity.attribute)
  for(j1 in 1:(length(ha)-1)){
    if( show.progress & ((j1 %% 10) == 0)) print(sprintf('completed:%.2f %%',j1*(j1-1)/2/total.iter*100))
    #system.time(
    for(j2 in (j1+1):length(ha)){
  	  ### get the mutual info between the 2 attributes
  	  e <- ha[j1]+ha[j2]-entropy(cbind(ent.att[,j1], ent.att[,j2]))
  	  ### or alternatively using the definition
  	  #p2 <- table((ent.att[,j1]*2+ ent.att[,j2]))
       #e2 <- ha[j1]+ha[j2]-sum(-p2/ne * log(p2/ne))
       mi[j1,j2] <- e
  	}
    #)	
  }
  if( show.progress ) print(sprintf('completed:%.2f %%',100))
}




      #write.table(round(mi,7),'mutual_info_GO_Yeast.xls',quote=FALSE,sep="\t")
      return(round(mi,7))
}
