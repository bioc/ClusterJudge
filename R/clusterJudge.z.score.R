clusterJudge.z.score <-
function(clusters, entity.attribute, nmb.randomizations=30){
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
     
#####  keep only attributes that have at least one gene
eat <- eat[,apply(eat,2,sum)>0]     
    if(dim(eat)[2] <=0 ) stop('No attributes found in the entity.attribute table for the clustered entities !')


##### z-score definition
#####   z = (MIreal − mean(MIrandom.unif.size))/srandom   (where s_random is the stdev of MI_random.unif.size)
#####   to maximize randomness take the same number of clusters but with aprox same size 
####  calculate z score and show plots 
no.clust      <- length(unique(clusters))  ### the number of clusters
RND.clust.size <- round(length(clusters)/no.clust) ### the size of the random clusters (as uniform as possible)
N.RND <- nmb.randomizations   ### number of randomizations
MIr <- vector(mode='numeric', length= nmb.randomizations)  ### will store the Mutual information after each randomization
clusters.RND <- rep(no.clust, length(clusters)) ###  create a randomized cluster of same size as the original
for(k in 1:length(clusters.RND))  clusters.RND[k] <- floor((k-1)/RND.clust.size) + 1  ### set the clusters in order because randomization will be later
#### add to last class until it's size is RND.clust.size-1
k=1; while(sum(clusters.RND == no.clust ) < (RND.clust.size-1) ){ clusters.RND[clusters.RND ==k][1] <- no.clust; k <- k+1}

HCr  <- ncol(eat)*entropy(clusters.RND)
HAr  <- sum(apply(eat,2,entropy))
cat('randomizing clusters ')
for(j in 1:N.RND){
     HCAr <- sum(apply(eat,2,function(v) entropy(cbind(sample(clusters.RND),v)))) ### here we randomize by sampling
     MIr[j]  <- HCr+HAr-HCAr
     cat('.') 	
}	

    calculate.mi <- function(cl,ea){
    	  n <- ncol(ea)
    	  n*entropy(cl)+sum(apply(ea,2,entropy)) - sum(sapply(1:ncol(ea), function(j) entropy(cbind(cl, ea[,j]))) )
    	}
    	
    MI <- calculate.mi(cl=clusters, ea=eat)

z.cl    <- (MI -mean(MIr))/ sd(MIr)
z.RND   <- (MIr -mean(MIr))/ sd(MIr)

 plot(bwplot( z.RND ~ rep('Random clusters of uniform size',length(MIr))
      ,ylim=c( min(c(z.cl,z.RND)), max(c(z.cl,z.RND)) )*1.1
      ,ylab='z-score = (MI.experiment - mean(MI.random))/stdev(MI.random)'
      ,panel=function(x,y,...){
      	panel.abline(v=1,h=seq(-10,100,10),alpha=0.5,lty=3)
      	panel.bwplot(x,y,...)
      	panel.xyplot(x,y,...,pch=19,alpha=0.2,jitter.x=TRUE)
      	panel.xyplot(1,z.cl,pch=19,alpha=1,col='red',cex=2)
      	panel.text(1,z.cl,pch=19,paste('z-score Clustering under study:',round(z.cl,2)), cex=0.8, pos=4)
      })
 )     

return(data.frame(z.score=c(z.cl, z.RND), method=c('experiment',rep('randomization',length(z.RND)))  ) )     
}
