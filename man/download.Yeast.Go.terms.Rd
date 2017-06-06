\name{download.Yeast.Go.terms}
\alias{download.Yeast.Go.terms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  downloads the description of Gene Ontology attributes
}
\description{
  The full description of Gene Ontology attributes are provided for example by the 
  Saccharomyces Genome Database (SGD).
}
\usage{
download.Yeast.Go.terms(url.GO.terms = "http://downloads.yeastgenome.org/curation/literature/go_terms.tab")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{url.GO.terms}{
     The web site address where the description of Gene Ontology attribute files is located.
}
}
\details{
  The file is expected to be in tab delimited format having the following columns: 'GOID' ,'GO_Term', 'GO_Aspect', 'GO_Term_Definition'
}
\value{
 a data.frame containing the information provided by the Gene Ontology consortium 
 
}
\references{
   http://www.yeastgenome.org/ http://www.geneontology.org/
}
\author{
  Adrian Pasculescu
}
\note{
  the function is optional. It only helps user in understanding what is the meaning of the GO attributes for the Yeast examples
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% download.Yeast.GO.mapping
\code{\link{download.Yeast.GO.mapping}}
}

\examples{

    GO.terms <- download.Yeast.Go.terms()
    str(GO.terms)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
