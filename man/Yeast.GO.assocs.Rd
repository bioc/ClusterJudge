\name{Yeast.GO.assocs}
\alias{Yeast.GO.assocs}
\docType{data}
\title{
    Gene Ontology attributes associated to Yeast Gene entities
}
\description{
    Gene Ontology attributes associated to Yeast Gene entities
}
\usage{data("Yeast.GO.assocs")}
\format{
    A data frame with 70487 observations on the following 2 variables.
    \describe{
    \item{\code{SGDID}}{a character vector}
    \item{\code{GOID}}{a character vector}
    }
}
\value{
    This data loads a  data frame with two columns: Yeast gene ids
    and their Gene Ontology ids
}
\references{
    Ashburner et al. Gene ontology: tool for the unification of biology 
    (2000) Nat Genet 25(1):25-9.

    The Gene Ontology Consortium. Gene Ontology Consortium: going forward. 
    (2015) Nucl Acids Res 43 Database issue D1049-D1056.
}
\examples{
data(Yeast.GO.assocs)
str(Yeast.GO.assocs) 
}
\keyword{datasets}
