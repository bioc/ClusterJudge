\name{validate.association}
\alias{validate.association}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 validates the associations between entities and attributes
}
\description{
  verifies if the input structure has two columns, if there are no NAs or NULLs
  and if thre are no duplicated associations
}
\usage{
validate.association(entity.attribute, message = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{entity.attribute}{
     data frame or matrix with 2 columns 
	The assumption is that first column represent some `entities` like gene names or gene ids. 
  And the second column represents `attributes` of entities (for example Gene Ontology ID 'GO:0007260' 
    which is `tyrosine phosphorylation of STAT protein`)

}
  \item{message}{
  logical. If set to TRUE will print a message in case of succesfull validation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Yeast.GO.assocs)

validate.association(Yeast.GO.assocs) 
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
