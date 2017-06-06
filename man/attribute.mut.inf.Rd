\name{attribute.mut.inf}
\alias{attribute.mut.inf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculates the mutual information between each attribute of the entity.attribute pairs
}
\description{
  calculates the mutual information based on the entropy. 
   The mutual information of the pair of attributes A,B is mi(A,b) = H(A)+H(B) -H(cbind(A,B))
   where H is the entropy.
}
\usage{
attribute.mut.inf(entity.attribute, show.progress = FALSE, alternative.calc = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{entity.attribute}{
     data frame or matrix with 2 columns 
	The assumption is that first column represent some `entities` like gene names or gene ids. 
  And the second column represents `attributes` of entities (for example Gene Ontology ID 'GO:0007260' 
    which is `tyrosine phosphorylation of STAT protein`) 
}
  \item{show.progress}{
     if set to TRUE will try to show periodically the estimated percentage of completion.
     Note that calculation of mutual information between Gene Ontology (GO) attributes
     takes a long time. We already have precalculated the GO mutual information
     for Yeast 
}
  \item{alternative.calc}{
    logical that if set to TRUE will try to use a direct calculation of the entropy from its definition
}
}
\details{
  
}
\value{
  a matrix having as rownames and columnnames the attribute names from the input entity.attribute structure.
  The values are the mutual informations between the attributes in the region above the main diagonal.  
  All the values below the diagonal are NA or zero.
  
}
\references{
     Gibbons, F.D. and Roth F.P., (2002)
    Judging the Quality of Gene Expression-Based Clustering Methods Using Gene Annotation. 
    Genome Research, vol. 12, pp1574-1581.

    Cover,T.M. and Thomas,J.A. 1991. Elements of information theory (ed. D.L. Schilling). Wiley-Interscience,New York.
}
\author{
  Adrian Pasculescu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mi.GO.Yeast}}
\code{\link{consolidate.entity.attribute}}
}
\examples{
data(mi.GO.Yeast) ### is a precalculated example of mutual information for the Yeast genes and their Gene Ontollogy (GO) attributes
str(mi.GO.Yeast)

#### the following example will generate the mutual information of 500 sampled pairs of Yeast , GO associations 
####    (otherwise it might take take minutes or hours to calculate the mutual information)
data(Yeast.GO.assocs) ### 
entity.attribute.sampled <- Yeast.GO.assocs[sample(1:nrow(Yeast.GO.assocs),500),]
mi.GO.Yeast.sampled <- attribute.mut.inf(entity.attribute.sampled, show.progress=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Yeast.GO.assocs }% use one of  RShowDoc("KEYWORDS")
\keyword{ mi.GO.Yeast }% __ONLY ONE__ keyword per line
\keyword{ validate.association }% __ONLY ONE__ keyword per line
