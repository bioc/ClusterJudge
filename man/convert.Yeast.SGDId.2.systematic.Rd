\name{convert.Yeast.SGDId.2.systematic}
\alias{convert.Yeast.SGDId.2.systematic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  converts the SGD Gene Id into the systematic name of the gene
}
\description{
  converts the Saccharomyces Genome Database (SGD) Gene Id 
  into the systematic name of the gene
}
\usage{
convert.Yeast.SGDId.2.systematic(SGD.Ids = c("S000007287", "S000007287", "S000007287"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SGD.Ids}{
%%     ~~Describe \code{SGD.Ids} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a matrix of characters containing the input SGD gene Ids as first column 
  and the corresponding Yeast systematic gene names as seccond column.
}
\references{
Yeast - Systematic Names Assignment in
http://www.yeastgenome.org/help/community

GF Berriz and FP Roth. The Synergizer service for translating gene, protein, 
and other biological identifiers. Bioinformatics 24(19)-2272 (2008).

}
\author{
  Adrian Pasculescu
}
\note{
  it requires online access to http://llama.mshri.on.ca/cgi/synergizer/serv
  where the conversion is performed.
  
  It uses the httr and jsonlite R packages.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

converted.genes <- convert.Yeast.SGDId.2.systematic(
                   SGD.Ids = c("S000007287", "S000007287", "S000007287"))

str(converted.genes)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
