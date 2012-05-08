\name{elastic}
\alias{elastic}
\title{
Elasticities of EASI demand system
}
\description{
This function calculates price elasticities, demographic elasticities and income elasticities.
}
\usage{
elastic(object = object, type = c("price", "income", "demographics"), sd = FALSE)
}
\arguments{
  \item{object}{
an object of type easi which contains all the information about the previous estimates. Here, in particular, the estimates of the easi system.
}
  \item{type}{
type can take values "price", "demographics" or "income" if one wants to calculate price elasticities, demographics elasticities or income elasticities, respectively.
}
  \item{sd}{
sd (for standard deviations) must be specified as TRUE if one wants to calculate the standard deviations (FALSE otherwise).
}
}
\details{
The equations for the budget shares elasticities are presented in Pendakur and Lewbel. The elasticities for quantities are derived in the context of the easi package.
}
\value{
 \item{EP}{semi-elasticities of budget shares with respect to prices.}
\item{ER}{semi-elasticities of budget shares with respect to real expenditures.}
\item{EZ}{semi-elasticities of budget shares with respect to demographics.}
\item{EPS}{matrix of compensated quantity derivatives with respect to unlogged prices.}
\item{EPQ}{compensated (good-specific) expenditures with respect to prices.}
 \item{EP_SE}{standard deviations of semi-elasticities of budget shares with respect to prices.}
\item{ER_SE}{standard deviations of semi-elasticities of budget shares with respect to real expenditures.}
\item{EZ_SE}{standard deviations of semi-elasticities of budget shares with respect to demographics.}
\item{ELASTPRICE}{elasticities of quantities with respect to prices.}
\item{ELASTINCOME}{elasticities of quantities with respect to income.}
\item{ELASTPRICE_SE}{standard deviations of Elasticities of quantities with respect to prices.}
\item{ELASTINCOME_SE}{standard deviations of Elasticities of quantities with respect to income.}
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

#******* price elasticities *******************
#******* without standard deviations **********
ep1 <- elastic(est,type="price",sd=FALSE)

#******* price elasticities *******************
#******* with standard deviations *************
ep2 <- elastic(est,type="price",sd=TRUE)

#******* demographic elasticities *************
#******* without standard deviations **********
ed <- elastic(est,type="demographics",sd=FALSE)

#******* Elasticities income ******************
#******* without standard deviations **********
ei <- elastic(est,type="income",sd=FALSE)
}
}
\keyword{EASI}
\keyword{Elasticities}
