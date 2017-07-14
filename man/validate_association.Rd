\name{validate_association}
\alias{validate_association}
\title{
    validates the associations between entities and attributes
}
\description{
    verifies if the input structure has two columns, if there are no NAs or NULLs
    and if thre are no duplicated associations
}
\usage{
validate_association(entity.attribute, message = TRUE)
}
\arguments{
    \item{entity.attribute}{
    data frame or matrix with 2 columns 
	  The assumption is that first column represent some `entities` like gene names or gene ids. 
    And the second column represents `attributes` of entities (for example 
    Gene Ontology ID 'GO:0007260' which is `tyrosine phosphorylation of STAT protein`)
}
    \item{message}{
    a logical. If set to TRUE will print a message in case of succesfull validation.
}
}
\value{
    Returns TRUE and print "Validation OK!" message in case of scuccess.
    Otherwise stops with an explanatory message.
}
\references{
    Gibbons, F.D. and Roth F.P., (2002)
    Judging the Quality of Gene Expression-Based Clustering Methods Using Gene Annotation. 
    Genome Research, vol. 12, pp1574-1581.}
\author{
    Adrian Pasculescu
}
\examples{
data(Yeast.GO.assocs)

validate_association(Yeast.GO.assocs) 
  
}
