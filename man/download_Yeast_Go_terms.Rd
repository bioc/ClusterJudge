\name{download_Yeast_Go_terms}
\alias{download_Yeast_Go_terms}
\title{
    downloads the description of Gene Ontology attributes
}
\description{
    The full description of Gene Ontology attributes are provided for example
    by the Saccharomyces Genome Database (SGD).
}
\usage{
download_Yeast_Go_terms(
  url.GO.terms = "http://downloads.yeastgenome.org/curation/literature/go_terms.tab")
}
\arguments{
    \item{url.GO.terms}{
    The web site address where the description of Gene Ontology attribute files is located.
}
}
\details{
    The file is expected to be in tab delimited format having the following columns: 
    'GOID' ,'GO_Term', 'GO_Aspect', 'GO_Term_Definition'
}
\value{
    a data.frame containing the information provided by the 
    Gene Ontology consortium 
}
\references{
    http://www.yeastgenome.org/ http://www.geneontology.org/
}
\author{
    Adrian Pasculescu
}
\note{
    the function is optional. It only helps user in understanding what is the 
    meaning of the GO attributes for the Yeast examples.
}
\seealso{
\code{\link{download_Yeast_GO_mapping}}
}

\examples{
    GO.terms <- download_Yeast_Go_terms()
    str(GO.terms)
}
