\name{simulations}
\alias{simulations}
\title{
Simulation tools in the context of the EASI demand system
}
\description{
The simulations function uses estimates of the easi function in order to simulate the new fitted budget shares and the new elasticities after changes in prices, demographics and total expenditure.
}
\usage{
simulations(object = object, log.price_new = log.price_new, var.soc_new = var.soc_new, 
             log.exp_new = log.exp_new)
}
\arguments{
  \item{object}{
an object of type easi which contains all the information about the previous estimation. Here, in particular, the estimates of the easi system.
}
  \item{log.price_new}{
the new price matrix (in logarithms).
}
  \item{var.soc_new}{
the new demographic variables matrix.
}
  \item{log.exp_new}{
the new total expenditures (in logarithms).
}
}
\details{
Simulation schemes are derived in the context of the easi package.
}
\value{
\item{summary}{Summary of estimation.}
\item{fitted.w}{Matrix of fitted budget shares.}
\item{CoefCov}{Covariance matrix of estimates.}
\item{y}{Implicit utility.}
  \item{shares}{
the matrix of budget shares.
}
  \item{log.price}{
the matrix of prices (in logarithms).
}
  \item{var.soc}{
the matrix of demographic variables.
}
  \item{log.exp}{
the logarithm of total expenditures.
}
  \item{y.power}{
the highest power of y (implicit utility) in the system.
}
  \item{neq}{
The number of goods (Excluding the last good).
}
  \item{nsoc}{
The number of demographic variables.
}
 \item{py.inter}{a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the price variables and y (implicit utility).
}
 \item{zy.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the demographic variables and y (implicit utility).
}
  \item{pz.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the prices and the demographic variables.
}
  \item{interpz}{
a vector which allows to choose the demographic variables to interact with the price. For example, interpz=c(3) means that price are interacted with the third demographic variable while interpz = c (1:n) means that prices are interacted with the first n demographic variables.
}
  \item{labels.shares}{a string which contains the names of budget shares. For example, labels.shares = c ("foods", "non foods", "others") indicates that the names of the three shares in the budget system are "foods", "non foods" and "others", respectively.
}
  \item{labels.price}{
a string which contains the names of log.prices. For example, labels.price = c ("pfoods", "pnonfoods", "pothers") indicates that the names of the three prices variables in the budget system are: "pfoods", "pnonfoods" and "pothers", respectively.
}
  \item{labels.soc}{
a string which contains the names of demographic variables. For example, labels.soc = c ("age", "gender", "hsize") indicates that the names of the three demographic variables in the budget system are "age", "non foods" and "hsize", respectively.
}
\item{interact}{a logical variable which is fixed to TRUE if at least one of the variables in either zy.inter, py.inter or pz.inter is fixed to TRUE.}
\item{fit3sls}{object of class systemfit.}
\item{varlist}{list of all variables.}
\item{dim_varlist}{number of variables.}

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

#******* Budget shares Matrix ****************
shares_HIX=hixdata[,2:10]
#****** price Matrix (in logarithms) *********
log.price_HIX=hixdata[,11:19] 
#****** demographic matrix *******************
var.soc_HIX=hixdata[,21:25]
#****** logarithm of total expenditure *******
#****** (here divised by a price index) ******
log.exp_HIX=hixdata[,20]  

#****** estimation ***************************
est <- easi(shares=shares_HIX,log.price=log.price_HIX,
             var.soc=var.soc_HIX,log.exp=log.exp_HIX)

#****** Changes in prices ********************
#****** PRIX.HIX --> PRIX.HIX.SIM1 = 0 *******
log.price_HIX.SIM1 <- log.price_HIX
for (i in 1:ncol(log.price_HIX))
log.price_HIX.SIM1[,i] <- 0

#****** Changes in demographic variables *****
#****** var.soc.HIX --> var.soc.HIX.SIM1 = 0 *
var.soc_HIX.SIM1 <- var.soc_HIX
for (i in 1:ncol(var.soc_HIX))
var.soc_HIX.SIM1[,i] <- 0

#******* simulation **************************
sim <- simulations(est,log.price_new=log.price_HIX.SIM1,
       var.soc_new=var.soc_HIX.SIM1,log.exp_new=log.exp_HIX)

#******* corresponding Engel curves **********
eng2 <- engel(sim,file="essai2",sd=FALSE)
}
}
\keyword{EASI}
\keyword{Simulation}
