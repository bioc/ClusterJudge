\name{mi.GO.Yeast}
\alias{mi.GO.Yeast}
\docType{data}
\title{
precalculated mutual information between Gene Ontology attributes of Yeast genes
}
\description{
precalculated mutual information between Gene Ontology attributes of Yeast genes
}
\usage{data("mi.GO.Yeast")}
\format{
    The format is:
    num [1:2266, 1:2266] NA NA NA NA NA NA NA NA NA NA ...
    - attr(*, "dimnames")=List of 2
     ..$ : chr [1:2266] "GO:0000001" "GO:0000002" "GO:0000009" "GO:0000011" ...
     ..$ : chr [1:2266] "GO:0000001" "GO:0000002" "GO:0000009" "GO:0000011" ...
}
\details{
    for convenienece this data set was pre-generated using attribute.mut.inf from this package
}
\value{
    This data loads a symetric matrix of mutual information values
    calculated between pairs of Gene Ontology attributes of Yeast genes
}
\references{
Wikipedia: Mutual Information https://en.wikipedia.org/wiki/Mutual_information
}
\seealso{
    \code{\link{attribute.mut.inf}}
}
\examples{
data(mi.GO.Yeast)
str(mi.GO.Yeast) 

}
\keyword{datasets}
