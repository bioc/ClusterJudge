entities.attribute.stats <-
function(entity.attribute
                                 ,  min.entities.per.attr=NULL
                                 , entity.space.name='Yeast genes'
                                 , attribute.space.name='Gene Ontology' ){
    entities.per.attr <- tapply(entity.attribute[,1],entity.attribute[,2],length)
    ### the entities.per.attr are distributed exponentially see table below
    entities.per.attr.tab <- table(entities.per.attr)
    if(is.null(min.entities.per.attr)){
      ### so we pick as threshold the values for wich the frequency is higher than 1/3 of the max value
      filt <- (entities.per.attr.tab <= max(entities.per.attr.tab, na.rm=TRUE)/3  )  
      min.entities.per.attr = min(as.numeric(names(entities.per.attr.tab)[filt]), na.rm=TRUE)  ### this is the proposed min entities per attribute 
    }
    if(!is.numeric(min.entities.per.attr)){
       stop(' min.entities.per.attr must be NULL or numeric !')
    }	    
    attr.ignored<- sum(entities.per.attr <= min.entities.per.attr)
    attr.kept   <- length(entities.per.attr)-attr.ignored
    pct.kept    <- attr.kept/length(entities.per.attr)*100
    plot(histogram(~entities.per.attr
           ,scales=list(x=list(log=10)), pch=1, cex=0.3, alpha=0.9
           ,par.settings=list(par.sub.text=list(cex=0.7))
           ,xscale.components = xscale.components.log10ticks
           ,xlab=paste('Number of',entity.space.name,'per', attribute.space.name,'attribute')
           ,main=paste('Distribution of', entity.space.name,'per', attribute.space.name,'attribute\n and the proposed fraction of attributes to be kept')
           ,panel=function(x,...){
           	  panel.histogram(x,...)
           	  panel.rect(-10,-10,log(min.entities.per.attr,10),100, border='transparent',col='lightgrey', alpha=0.4)
           	  panel.abline(v=log(min.entities.per.attr,10),col='red')
           	  panel.text(log(min.entities.per.attr,10),0,paste('# Attributes ignored:\n',attr.ignored),cex=0.8,pos=2)
           	  panel.text(log(min.entities.per.attr,10),0,paste('# Attributes kept:\n',length(entities.per.attr)-attr.ignored, '(',round(pct.kept,2),'%)')
           	           ,cex=0.8, pos=4)
           	}
           ,sub=paste('Proposed threshold for attributes present in less than:', min.entities.per.attr,'entities')
           )	
         )
    return( min.entities.per.attr)           	
}
