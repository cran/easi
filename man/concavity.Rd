\name{concavity}
\alias{concavity}
\title{
Concavity of EASI cost function
}
\description{
This function provide a test of local concavity
}
\usage{
concavity(object = object)
}
\arguments{
  \item{object}{
an object of type easi which contains all the information about the previous estimates. Here, in particular, the estimates of the easi system.
}
}
\details{
A semi-negative definite Hessian matrix is a necessary and sufficient condition to consider that the cost function is concave. The calculation of Hessian matrix and the test of local concavity is performed in the context of the easi package.
}
\value{
\item{a}{A string that indicates the result of the test}
}
\references{
Arthur Lewbel and Krishna Pendakur "Tricks with Hicks : The EASI Demand System"
The American Economic Review 99(3), 827-863, 2009.

Krishna Pendakur "EASI made Easier",
www.sfu.ca/pendakur/EASI made Easier.pdf, 2008.
}
%\author{}
%\note{}
\examples{
\donttest{
data(hixdata)

#******* Budget shares Matrix *****************
shares_HIX=hixdata[,2:10]
#******* price Matrix (in logarithms) *********
log.price_HIX=hixdata[,11:19] 
#******* demographic matrix *******************
var.soc_HIX=hixdata[,21:25]
#******* logarithm of total expenditure *******
#******* (here divised by a price index) ******
log.exp_HIX=hixdata[,20]  

#******* estimation ***************************
est <- easi(shares=shares_HIX,log.price=log.price_HIX,
            var.soc=var.soc_HIX,log.exp=log.exp_HIX)

#******* test of local concavity **************
concavity(est)}

}
\keyword{EASI}
\keyword{Elasticities}
