clusterJudge <-
function(clusters,entity.attribute, plot.notes=''){
if(is.null(names(clusters))) stop('The clusters must have as names the entity names of the reference entity.attribute table (e.g. Yeast gene names)')
#####  calculate MI between Clusters and each Attribute: MI(C,Ai) . 
#####     For all attributes (assuming independence) caculate the sum MI(C,A) = Sum(MI(C,Ai)) = n*E(C)+sum(E(Ai))-sum(E(C,Ai))
#####         where C is the variable that has the cluster id ! 
#####  we will restrict the GO attributes table to the entities present in the clusters
    if(!validate.association(entity.attribute, message=FALSE)) return() 

    ##### drop clustered entities not in entities and reciprocally
    common.entities <- intersect(names(clusters), entity.attribute[,1])
    clusters         <- clusters[common.entities]
    entity.attribute <- entity.attribute[entity.attribute[,1] %in% common.entities,]
     
    eat <- table(entity.attribute)
      
    calculate.mi <- function(cl,ea){
    	  n <- ncol(ea)
    	  n*entropy(cl)+sum(apply(ea,2,entropy)) - sum(sapply(1:ncol(ea), function(j) entropy(cbind(cl, ea[,j]))) )
    	}
    	
    MI <- calculate.mi(cl=clusters, ea= eat)

    #### completely randomize
    MI.RND.TOTAL <- calculate.mi(cl=sample(clusters), ea=eat)

#### select some 50 levels of randomization
    # use power of 2 as number of swaps
    n.rnd <- 12
    mi.RND <- vector(mode='numeric', length=n.rnd)
    mi.ea <- sum(apply(eat,2,entropy))
    n <- ncol(eat) 
    cl <- clusters
    cl.ids <- names(table(cl))
    ### the entropy of the swapped clusters does NOT change so we cakculate it just once
    mi.cl <- n*entropy(cl)
    cat('randomizing clusters ')
    for(r in 1:n.rnd){
    	for(i in 1:(2^(r-1)) ){   
   		c1.2 <- sample(cl.ids,2)	
   		g1 <- sample(names(cl)[cl==c1.2[1]],1)
   		g2 <- sample(names(cl)[cl==c1.2[2]],1)
   		cl[g1] <- c1.2[2]
   		cl[g2] <- c1.2[1]
    	}	
        cat(".")
    	mi.RND[r] <-  mi.cl + mi.ea -  sum(sapply(1:ncol(eat), function(j) entropy(cbind(cl, eat[,j]))) )
    }	
    cat("\n")
    
    
    #### plot result
    plot(xyplot(c(MI,mi.RND[1:n.rnd],MI.RND.TOTAL)/MI~c(0:n.rnd,n.rnd+2)
      ,grid=TRUE,  cex=1.6, type=c('b')
      ,xlab='Number of random swaps between original clusters' 
      ,ylab='Mutual Information - as a fraction of initial value'
  ,scales=list(x=list(cex=0.7, rot=45, at=c(0:(n.rnd+4)),labels=c(0,2^c(0:(n.rnd-1)),'...','full\nshuffle')))
      ,pch=19
      ,panel=function(x,y,...){
      	panel.xyplot(x,y,...)
      }	
      ,main='Decrease of relative Mutual Information\n after random swaping between pairs of clusters'
      ,sub=plot.notes
      )
    )

   data.frame(mut.inf=c(MI,mi.RND[1:n.rnd] ,MI.RND.TOTAL)
             ,nmb.of.swaps = c(0,2^(0:(n.rnd-1)),' full shuffle'))
}
