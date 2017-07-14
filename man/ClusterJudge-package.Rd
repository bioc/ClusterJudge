\name{ClusterJudge-package}
\alias{ClusterJudge-package}
\alias{ClusterJudge}
\docType{package}
\title{
\packageTitle{ClusterJudge}
}
\description{
\packageDescription{ClusterJudge}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{ClusterJudge}
\packageIndices{ClusterJudge}
    ClusterJudge is a way to judge the quality of clustering performed elsewhere on some entities.
    This judgement is based on some additional entitity-attribute information.
    For example, using Gene Ontology annotated attributes offered by Saccharomyces Genome Database (SGD),
    it can judge the clusters of Yeast genes (i.e. entities) resulted from experiments related to the mitotic Yeast cell cycle.   
    This is done by evaluating the mutual information between a gene membership in a cluster, and the attributes it possesses. 
}
\author{
\packageAuthor{ClusterJudge}

Maintainer: \packageMaintainer{ClusterJudge}
}
\references{
    
    Gibbons, F.D. and Roth F.P., (2002)
    Judging the Quality of Gene Expression-Based Clustering Methods Using Gene Annotation. 
    Genome Research, vol. 12, pp1574-1581.
    
    Tavazoie, et al. (1999) Systematic determination of genetic network architecture. Nature Genetics, 22, 281-285
    
    
    Cho, et al. (1998) A genome-wide transcriptional analysis of the mitotic cell cycle. Molecular Cell,2, 65-73.
    
    Cover,T.M. and Thomas,J.A. 1991. Elements of information theory (ed. D.L. Schilling). Wiley-Interscience,New York.
    
}
\keyword{ package }
