\name{download.Yeast.GO.mapping}
\alias{download.Yeast.GO.mapping}
\title{
    downloads the Gene Ontology attributes mapped to Yeast genes
}
\description{
    The Gene Ontology attributes are provided for example by the Saccharomyces 
    Genome Database (SGD) as a zipped file with a wealth of information. 
    The function filters from the file only the minimal information
    necessary: THe Gene Ontology attribute id (GOID) and the Yeast Gene Id 
    (SGDID) for which the attribute applies   
}
\usage{
download.Yeast.GO.mapping(
    yeast.GO.url = "http://downloads.yeastgenome.org/curation/literature/gene_association.sgd.gz")
}
\arguments{
    \item{yeast.GO.url}{
    The web site address where the Yeast genes and their Gene Ontology 
    attribute files is located.
}
}
\details{
    Only the unique associations between genes and their attributes are kept. 
    In the original downloaded file there might be redundanies due to the different 
    types of evidences used in the associations.

    The downloaded file is expected to be in a tab delimited format with comment lines 
    that start with the exclamation (!) caharacter.
    The commented lines are ignored. The expected column names can be obrtained
    by looking at the function code.  
}
\value{
    a data.frame with two columns:
    SGDID - the ids of the Yeast gene names (these are the `entities`)
    GOID  - the ids of the Gene Ontology attributes corresponding to the above genes
}
\references{
http://www.yeastgenome.org/ http://www.geneontology.org/
}
\author{
    Adrian Pasculescu
}
\note{
    since the download is time consuming, this package provides also as a dataset 
    called Yeast.GO.assocs containing the associations already downloaded 
    as of on Fenruary 2017.
}
\seealso{
\code{\link{Yeast.GO.assocs}} 
\code{\link{help}}, ~~~
}
\examples{

Yeast.GO.assocs <- download.Yeast.GO.mapping()

}
