\name{easi}
\alias{easi}
\alias{coef.easi}
\alias{residuals.easi}
\alias{summary.easi}
\alias{predict.easi}
\alias{vcov.easi}
\title{
Exact Affine Stone Index (EASI) Demand System Estimation
}
\description{
This function estimates the Exact Affine Stone Index (EASI) Demand System of Pendakur and Lewbel (2009).
}
\usage{
easi(shares = shares, log.price = log.price, var.soc = var.soc, 
     log.exp = log.exp, y.power=FALSE, labels.share=FALSE, labels.soc=FALSE, 
     py.inter=FALSE, zy.inter=FALSE, pz.inter=FALSE, interpz=FALSE)
\method{coef}{easi}(object,...)
\method{residuals}{easi}(object,...)
\method{summary}{easi}(object,...)
\method{predict}{easi}(object,...)
\method{vcov}{easi}(object,...)
}
\arguments{
  \item{shares}{
the matrix of budget shares.
}
  \item{log.price}{
the matrix of log.prices (in logarithms).
}
  \item{var.soc}{
the matrix of demographic variables.
}
  \item{log.exp}{
the logarithm of total expenditure.
}
  \item{y.power}{
the highest power of y (implicit utility). Default is y.power=3.
}
  \item{labels.share}{a string which contains the names of budget shares. For example, labels.share = c ("foods", "non foods", "others") indicates that the respective names of the three shares in the budget system are "foods", "non foods" and "others".
}
  \item{labels.soc}{
a string which contains the names of demographic variables. For example, labels.soc = c ("age", "gender", "hsize") indicates that the respective names of the  three demographic variables in the budget system are "age", "non foods" and "hsize".
}
  \item{py.inter}{a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the price variables and y (implicit utility).
}
  \item{zy.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the demographic variables and y (implicit utility).
}
  \item{pz.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between prices and demographic variables.
}
  \item{interpz}{
a vector which allows to choose which demographic variables to interact with prices. For example, interpz=c(3) means that prices are interacted with the third demographic variable while interpz = c (1:n) means that prices are interacted with the first n demographic variables.
}
  \item{object}{
an object of type easi which contains all informations about the previous estimates. Here, in particular, the estimates of the easi system.
}
  \item{...}{
others options.
}
}
\details{
The EASI demand system, in the form of budget shares, is given by (for each good \eqn{j \in [1:J]}):
\cr

\eqn{w^j = \sum_{r=1}^R b_r^j y^r  + \sum_{t=1}^T g_t^j z_t + \sum_{k=1}^T \sum_{t=1}^T a_{jkt} z_t \ln p^k  + \sum_{k=1}^J b_{jk} \ln p^k y + \sum_{t=2}^T h_{t}^j z_t y + \varepsilon_j},
\cr

where \eqn{w} is the matrix of budget shares, \eqn{R} is the hightest power of implicit utility \eqn{y}, \eqn{T} denotes the number of demographic variables, \eqn{z} is the matrix of demographic variables, \eqn{J} is the number of goods and \eqn{p} is the matrix of prices. Implicit utility \eqn{y} is given by:
\cr

\eqn{ \phantom{123456789} {\displaystyle y = \frac{\ln x - \sum_{j=1}^J w_j \ln p^j+ 1/2 \sum_{j=1}^J \sum_{k=1}^J a_{jkt} z_t \ln p^j \ln p^k}{1-\frac{1}{2}\sum_{j=1}^J \sum_{k=1}^J b_{jk} \ln p^j \ln p^k}}\\[16pt].}
\cr

The system is estimated by an iterated linear 3SLS. This estimator is similar to the estimator suggested by Blundell and Robin (1999) for the QUAIDS, except for the treatment of the endogeneity of \eqn{w_j \ln p_j}.
}
\value{
  

  an object of class \code{c("easi", "systemfit")}.
    
  A \code{"easi"} object has the following elements :
  
\item{summary}{summary of estimation.}
\item{fitted.w}{the matrix of fitted budget shares.}
\item{Residuals}{the matrix of residuals.}
\item{CoefCov}{the covariance matrix of the estimates.}
\item{y}{the implicit utility.}
\item{varlist}{the list of all variables.}
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
the logarithm of total expenditure.
}
  \item{y.power}{
the highest power of y (implicit utility) in the system.
}
  \item{neq}{
the number of goods (Excluding the last good).
}
  \item{nsoc}{
the number of demographic variables.
}
  \item{py.inter}{a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the price variables and y (implicit utility).
}
  \item{zy.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between the demographic variables and y (implicit utility).
}
  \item{pz.inter}{
a logical variable which is fixed to TRUE (FALSE otherwise) if the user wants to enable the interaction between prices and demographic variables.
}
  \item{interpz}{
a vector which allows to choose the demographic variables to interact with prices. For example, interpz=c(3) means that prices are interacted with the third demographic variable, while interpz = c (1:n) means that prices are interacted with the first n demographic variables.
}
 \item{labels.shares}{a string which contains the names of budget shares. For example, labels.shares = c ("foods", "non foods", "others") indicates that the names of the three shares in the budget system are "foods", "non foods" and "others", respectively.
}
  \item{labels.price}{
a string which contains the names of prices. For example, labels.price = c ("pfoods", "pnonfoods", "pothers") indicates that the names of the three log.prices variables in the budget system are: "pfoods", "pnonfoods" and "pothers", respectively.
}
  \item{labels.soc}{
a string which contains the names of the demographic variables. For example, labels.soc = c ("age", "gender", "hsize") indicates that the names of the  three demographic variables in the budget system are "age", "non foods" and "hsize", respectively.
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

#******* Budget shares Matrix ************
shares_HIX=hixdata[,2:10]
#******* price Matrix (in logarithms) ****
log.price_HIX=hixdata[,11:19] 
#******* Demographic matrix **************
var.soc_HIX=hixdata[,21:25]
#******* logarithm of total expenditure **
#******* (here divised by a price index) *
log.exp_HIX=hixdata[,20]  

#******* Estimation... *******************
est <- easi(shares=shares_HIX,log.price=log.price_HIX,
             var.soc=var.soc_HIX,log.exp=log.exp_HIX)

#******* Coefficients ********************
coef(est)

#******* Summary *************************
summary(est)

#******* fitted values *******************
head(predict(est))

#******* Residuals ***********************
head(residuals(est))

#******* Covariance Matrix ***************
dim(vcov(est))
}
}
\keyword{EASI}


