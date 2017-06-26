\name{convert.Yeast.SGDId.2.systematic}
\alias{convert.Yeast.SGDId.2.systematic}
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
\arguments{
    \item{SGD.Ids}{
  	a vector of characters. Each element is a string representing a gene.
}
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
\examples{

converted.genes <- convert.Yeast.SGDId.2.systematic(
                   SGD.Ids = c("S000007287", "S000007287", "S000007287"))

str(converted.genes)

}
