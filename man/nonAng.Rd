\name{nonAng}
\alias{nonAng}
\docType{data}
\title{
Non-Angiogenic Ovarian Cancer Network Edge List
}
\description{
PANDA network was infered by integrating non-angiogenic ovarian cancer microarray expression data with transcription factor binding sites. The algorithm, expression data, and the network prior was downloaded from https://sites.google.com/a/channing.harvard.edu/kimberlyglass/tools/panda.
}
\usage{data("nonAng")}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(nonAng)
head(nonAng)
input=cbind(nonAng,ang[,3])
}
\keyword{datasets}
