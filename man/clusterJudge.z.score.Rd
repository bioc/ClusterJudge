\name{clusterJudge.z.score}
\alias{clusterJudge.z.score}
\title{
    calculates the `cluster judge` z-score  
}
\description{
    calculates the `cluster judge` z-score as defined in the reference
    The z-score is based on shuffling the clusters at random and calculating the total mutual information relative to the entity.attribute table.
    After the selecetd number of randomizations the mean MR and standard deviation SDR of the mutual information is used in the
    definition of the z.score = (MI- MIR)/SDR  where MI is the mutual information of the original clustering. 
    The higher the z.score the better the clustering. 
    A box-and-wisker  plot is generated that shows how far is the clustering versus random clustering 
    based on the mutual information to the selected entitity.attribute 
}
\usage{
clusterJudge.z.score(clusters, entity.attribute, nmb.randomizations = 30)
}
\arguments{
    \item{clusters}{
    a named vectors of integers (or a factor). 
    The names (or the levels of the factor) must match some (as many as possible) of the rownames of the entity.attribute table.
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

    \item{nmb.randomizations}{
    number of randomization iterations
    }
}
\value{
    a data.frame
    with the number of randomization shuffles
    and the total mutual information calculated after each of the shuffles
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
    a dot is printed on the console after each randomization (shuffling) step
}
\examples{

library('yeastExpData')
data(ccyclered)

clusters <- ccyclered$Cluster
###  convert from Gene names to the new standard of Saccharomyces Genome Database (SGD) gene ids
ccyclered$SGDID <- sub('^S','S00',ccyclered$SGDID)
names(clusters) <- ccyclered$SGDID

data(Yeast.GO.assocs)  #### obtain associations and consolidate them at uncertainty level 0.001
Yeast.GO.assocs.cons <- consolidate.entity.attribute(entity.attribute = Yeast.GO.assocs
                                                   , min.entities.per.attr =3
                                                   , mut.inf=mi.GO.Yeast
                                                   , U.limit = c(0.001)) 

#### calculate z.scores for the associations consolidated at 0.001 Uncertainty level
z.scores <- clusterJudge.z.score(clusters
                        , entity.attribute = Yeast.GO.assocs.cons[[0.001]]
                        , nmb.randomizations=30)
}
\keyword{ clusterJudge }
\keyword{ Yeast.GO.assocs }
