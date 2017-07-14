\name{consolidate_entity_attribute}
\alias{consolidate_entity_attribute}
\title{
    removes the redundant attributes based on the mutual information between attributes
}
\description{
    Calculates the mutual information between every pair of attributes and the uncertainty level 
    (as defined in the reference ...)
    Plots the distribution of the uncertainty and the selected levels of uncertainty 
    that are used as filteres.
}
\usage{
consolidate_entity_attribute(entity.attribute, min.entities.per.attr, mut.inf=NULL
                           , U.limit = c(0.8, 0.6, 0.4, 0.2, 0.1, 0.01, 0.001)
                           , plot.saveRDS.file=NULL)
}
\arguments{
    \item{entity.attribute}{
    data frame or matrix with 2 columns 
	The assumption is that first column represent some `entities` like gene names or gene ids. 
    And the second column represents `attributes` of entities (for example Gene Ontology ID 'GO:0007260' 
    which is `tyrosine phosphorylation of STAT protein`)
    Usually this is a consolidated entity.attribute where the attributes with very low number of entities 
    or with very low mutual information have been removed (see consolidate_entity_attribute 
    and the definition of Uncertainty on attributes mutual information)   
} 
 \item{min.entities.per.attr}{
    a number or NULL     
}
    \item{mut.inf}{
    the mutual information square matrix before applying any filteres based on the uncertainty
    or NULL
}
    \item{U.limit}{
    a numerical vector 
}
    \item{plot.saveRDS.file}{
    if not NULL must be a string represented a file location where the plot will be saved as an RDS object.
    The plot can be then retrieved at any time using readRDS function.
}

}
\value{
    a data frame (if mut.inf argument is null) or a list of data frames.
    The data frames are the `consolidated` entity.attribute structures.
    Where the consolidation of attributes is based on a minimum number of entities per attribute
    or on the mutual information between attributes.  
}
\references{
    Gibbons, F.D. and Roth F.P., (2002)
    Judging the Quality of Gene Expression-Based Clustering Methods Using Gene Annotation. 
    Genome Research, vol. 12, pp1574-1581.
}
\author{
    Adrian Pasculescu
}
\examples{

data(Yeast.GO.assocs)
data(mi.GO.Yeast)

Yeast.GO.assocs.cons <- consolidate_entity_attribute(entity.attribute = Yeast.GO.assocs
                                                   , min.entities.per.attr =3
                                                   , mut.inf=mi.GO.Yeast
                                                   ,U.limit = c(0.01, 0.001)
                                                   ,plot.saveRDS.file='this_plot.rds') ### save also the plot

str(Yeast.GO.assocs.cons)
pdf('my_plot.pdf') ### place the plot saved by consolidate_entity_attribute into a pdf file
plot(readRDS('this_plot.rds')) ### 
dev.off()
}
\keyword{ attribute_mut_inf }
