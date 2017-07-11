consolidate_entity_attribute <-
function(entity.attribute, min.entities.per.attr, mut.inf=NULL,  U.limit = c(0.8, 0.6, 0.4, 0.2, 0.1, 0.01, 0.001), plot.saveRDS.file=NULL){	
   if(!validate_association(entity.attribute, message=FALSE)) return()
   if(!is.null(min.entities.per.attr)){
   	   if(!is.numeric(min.entities.per.attr)) stop('min.entities.per.attr must be a numeric value!')
   entities.per.attr <- tapply(entity.attribute[,1],entity.attribute[,2],length)
   if( (min.entities.per.attr<= 0) | (min.entities.per.attr>= max(entities.per.attr))) stop('min.entities.per.attr must be a positive number less than max(entities per attr)')
   entity.attribute <-entity.attribute[entity.attribute[,2] %in% names(entities.per.attr)[entities.per.attr > min.entities.per.attr],]
   }
   
   if(is.logical(mut.inf)){
   	  ## we skip applying the mutual information filter
   	  if(! mut.inf) return( entity.attribute)
   }	
   
   if(is.logical(mut.inf) | is.null(mut.inf)){
   	 warning('Will calculate first the mutual information now! Takes time!')
   	 mut.inf <- attribute_mut_inf(entity.attribute, show.progress=TRUE)
   }	
      
   if(length(intersect(colnames(mut.inf),entity.attribute[,2])) != ncol(mut.inf)) stop('Mutual Information and entity.attribute do not have same attribute names!')
   if(ncol(mut.inf) != nrow(mut.inf)) stop('Mutual information must be a square (upper triangular) matrix !')
   	
   mi.max <- max(as.vector(mut.inf),na.rm=TRUE)
   U <- mut.inf/mi.max
   U.names <- kronecker(rownames(U),t(colnames(U)),paste,sep=';')
    
   #### Show the histogram of U
   Uv <- as.vector(U)    ; U.names.v <- as.vector(U.names)
   filter.U <- !is.na(Uv)
   Uv <- Uv[filter.U]  ; U.names.v <- U.names.v[filter.U]
   #### set to same very small value all very small values
   Uv[Uv <= 3e-06] <- 3e-06
   #### get the max of histogram in percentages
   h <- hist(log(Uv,10),plot=FALSE)
   h.max.pct <- 0.6*max(h$counts/sum(h$counts,na.rm=TRUE))*100
   p1<-histogram(Uv, type='percent'	           
            ,scales=list(x=list(log=10),y=list(tick.number=11)), alpha=0.9
        ,xscale.components = xscale.components.log10ticks
        ,xlab='Uncertainty (MutInf/max(MutInf)\n(red: selected limits and number of removed elements)'
        ,ylab='Percent of Total (Number of attribute pairs) '
        ,panel=function(x,...){
        	panel.abline(h=seq(0,100,5),lty=3,alpha=0.8)
        	panel.histogram(x,...)
        	panel.abline(v=log(U.limit,10),col='red',alpha=0.8)
        	panel.rect(log(U.limit,10),-10,10,200, col='lightgrey',alpha=0.1)
        	panel.text(x=log(U.limit,10),y=seq(0,h.max.pct,length=length(U.limit))
        	          ,paste(sapply(U.limit,function(l) sum(Uv>=l) ),' pairs\n'
        	                ,'(about ',sapply(U.limit,function(l) round(length(unique(as.vector(sapply(U.names.v[Uv>=l],function(p) unlist(strsplit(p,';'))))))/2) ),' attribs)' )
        	          ,col='red',cex=0.6,pos=4, alpha=0.8)
        	panel.arrows(log(U.limit,10)    , seq(0,h.max.pct,length=length(U.limit))
        	            ,log(U.limit*4,10)  , seq(0,h.max.pct,length=length(U.limit)), col='red',alpha=0.6, length=0.07, lwd=2)         
        }	
        ,main=paste('Distribution of Uncertainty\n between',nrow(U)*(nrow(U)-1)/2,'pairs of',nrow(U),'attributes')
        )
   plot(p1)
   if(!is.null(plot.saveRDS.file)) tryCatch( {saveRDS(p1,plot.saveRDS.file)}, error=function(e){ print(e)})

   entity.attribute.list<-list()
   for(u.crt in U.limit){
      crt.pairs <- U.names.v[Uv>=u.crt]
      entity.attribute.crt <- entity.attribute
      ###  see in how many pairs is each attribute
      i <- 0
      while( length(crt.pairs)>0 ){
        crt.attribs <- unique(unlist(strsplit(crt.pairs,';')))
        crt.attribs.count <- table(crt.attribs)[order(table(crt.attribs),decreasing=TRUE)]
      	attr.1  <- names(crt.attribs.count)[1]
      	attr.2  <- sub(';','',sub(attr.1,'',grep(attr.1, crt.pairs,value=TRUE))) ### list of the second attributes pairing attr.1
      	### get rid of the second GO attribute from gene.attribute.u and from the crt.pairs
      	entity.attribute.crt <- entity.attribute.crt[entity.attribute.crt[,2] %in% setdiff(entity.attribute.crt[,2], attr.2), ]
     	crt.pairs <- grep(paste(attr.2,collapse='|'),crt.pairs,invert=TRUE,value=TRUE) ### ignore the pairs that have and attribute in the attr.2 array
      	if((i %% 10) == 0) cat('.'); i <- i+1
      }
      cat("\n");
      print(paste('For the uncertainty limit',u.crt,'removed:',length(unique(entity.attribute[,2]))-length(unique(entity.attribute.crt[,2])),'attributes'))
      entity.attribute.list[[as.character(u.crt)]]<-entity.attribute.crt
   }

   return(entity.attribute.list)	
}
