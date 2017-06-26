\name{clusterJudge}
\alias{clusterJudge}
\title{
    judges clustering using an entity.attribute table 
}
\description{
    calculates the sum of mutual information between clusters and each of the attributes
    For example if we note mutual information between Clusters and each Attribute i as: 
        MI(C,Ai)  (where C is the variable that has the cluster ids)
    then, for all attributes (assuming independence) it records the sum:
        MI(C,A) = Sum(MI(C,Ai)) = n*E(C)+sum(E(Ai))-sum(E(C,Ai))

    Then swaps elements between a random pair of clusters and repeats the above calculations.
    The swapping is repeated  2 times, 4 times, 8 times ... up to 2^12 times
    Finaly a complete shuffle of clustering is applied.
    The resulting variation of the total mutual information as a fraction of the initial value
    is ploted against the increasing number of swaps (and the final suffling)
    A god clustering should have a pronounced decrease of the total mutual information 
    when the of number of swaps is increased.
}
\usage{
    clusterJudge(clusters, entity.attribute, plot.notes = "")
}
\arguments{
    \item{clusters}{
    a named vectors of integers (or a factor). 
    The names (or the levels of the factor) must match some (as many as possible) 
    of the entities of the entity.attribute structure.
}
    \item{entity.attribute}{
    data frame or matrix with 2 columns 
	The assumption is that first column represent some `entities` like gene names or gene ids. 
    And the second column represents `attributes` of entities (for example Gene Ontology ID 'GO:0007260' 
    which is `tyrosine phosphorylation of STAT protein`)
    Usually this is a consolidated entity.attribute where the attributes with very low number of entities 
    or with very low mutual information have been removed (see consolidate.entity.attribute 
    and the definition of Uncertainty on attributes mutual information)   
}
    \item{plot.notes}{
    a string that will be added to the plot as explanation of what clustering represents.
}
}
\value{
    a data.frame
    with the number of swapps between clusters used for randomization
    and the total mutual information calculated after each of the sets of swaps
    The last value is for the full shuffle of the clusters. 
    Since the shuffle is using the base \code{sample} function witout setting a random seed, the last value will vary.
}
\references{
    Gibbons, F.D. and Roth F.P., (2002)
    Judging the Quality of Gene Expression-Based Clustering Methods Using Gene Annotation. 
    Genome Research, vol. 12, pp1574-1581.
}
\author{
    Adrian Pasculescu
}
\note{
    a dot is printed on the console after each of the 12 sets of swaps
}

\seealso{
    \code{\link{consolidate.entity.attribute}}
    \code{\link{help}}
}
\examples{

library('yeastExpData')
data(ccyclered)

clusters <- ccyclered$Cluster
###  convert from Gene names to the new standard of Saccharomyces Genome Database (SGD) gene ids
ccyclered$SGDID <- sub('^S','S00',ccyclered$SGDID)
names(clusters) <- ccyclered$SGDID
str(clusters)

data(Yeast.GO.assocs)
str(Yeast.GO.assocs)

### get the pre-calculated mutual information between all pairs of Gene Ontology (GO) attributes of Yeast genes
data(mi.GO.Yeast)  
str(mi.GO.Yeast)


### consolidate (reduce) the number of attributes by removing attributes with very few entities 
###    and removing ones very corelated (i.e. with high mutual information)
Yeast.GO.assocs.cons <- consolidate.entity.attribute(entity.attribute = Yeast.GO.assocs
                                                   , min.entities.per.attr =3
                                                   , mut.inf=mi.GO.Yeast
                                                   , U.limit = c(0.8, 0.6, 0.4, 0.2, 0.1, 0.001)) 
                                                   
#### apply clusterJudge for enntity attributes having a very low uncertainty (0.001)
mi.by.swaps<-clusterJudge(clusters, entity.attribute=Yeast.GO.assocs.cons[["0.001"]]
                         , plot.notes='Yeast clusters judged at uncertainty level 0.001 - Ref: Tavazoie S,& all 
`Systematic determination of genetic network architecture. Nat Genet. 1999`')  

}
\keyword{ attribute.mut.inf }
\keyword{ consolidate.entity.attribute }
