\name{equiv.income}
\alias{equiv.income}
\title{
Equivalent income of EASI demand system
}
\description{
The function estimates the equivalent income under the EASI system. Two periods are considered : the reference period (before changes in prices and/or total expenditures) and the current period (after changes in prices and/or total expenditures).

The function calculates the equivalent income, namely the necessary income in the reference period (with reference prices) in order to have the same implicit utility  in the current period (with current prices and current total expenditures).

}
\usage{
equiv.income(object = object, log.exp_ref = log.exp_ref, log.exp_cur = log.exp_cur, 
                 log.price_ref = log.price_ref, log.price_cur = log.price_cur)
}
\arguments{
  \item{object}{
an object of type easi which contains all the information about the previous estimation. Here, in particular, the estimates of the easi system.
}
  \item{log.exp_ref}{
the logarithm of the total expenditures of the reference period.
}
  \item{log.exp_cur}{
the logarithm of the total expenditures of the current period.
}
  \item{log.price_ref}{
the logarithm of the prices of the reference period.
}
  \item{log.price_cur}{
the logarithm of the prices of the current period.
}
}
\details{
Equivalent income is derived specially in the context of easi package.
}
\value{
 \item{equiv_income}{Equivalent Income.}
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
# We propose to simulate a "hybrid model" to 
#          illustrate how to use the function:

# This model is composed of five budget shares whose means are 
# respectively equal to 0.25, 0.15, 0.20, 0.30 and 0.10.
w1 <- rnorm(3000,mean=0.25,sd=0.05)
w2 <- rnorm(3000,mean=0.15,sd=0.05)
w3 <- rnorm(3000,mean=0.20,sd=0.05)
w4 <- rnorm(3000,mean=0.30,sd=0.05)
w5 <- 1-w1-w2-w3-w4 
shares_SIM=data.frame(w1,w2,w3,w4,w5)

# We simulate five price vectors, whose means are respectively 
#          equal to 25, 15, 20, 30 and 10:

p1 <- log(rnorm(3000,mean=25,sd=3))
p2 <- log(rnorm(3000,mean=15,sd=2))
p3 <- log(rnorm(3000,mean=20,sd=3))
p4 <- log(rnorm(3000,mean=30,sd=4))
p5 <- log(rnorm(3000,mean=10,sd=1))
log.price_SIM=data.frame(p1,p2,p3,p4,p5)

# We simulate four demographics variable : V1, V3, V4 who 
#   are dummy variables, and V2 who take its values in N+

V1 <- abs(round(rnorm(3000,mean=0.7,sd=0.2)))
V2 <- abs(round(rnorm(3000,mean=2,sd=1)))+1
V3 <- abs(round(rnorm(3000,mean=0.7,sd=0.2)))
V4 <- abs(round(rnorm(3000,mean=0.7,sd=0.2)))
var.soc_SIM=data.frame(V1,V2,V3,V4)


# Finally, we simulate a vector of total expenditure 
#         whose the average is 1200.

log.exp_SIM <- log(rnorm(3000,mean=1200,sd=200))

# Estimation 
est <- easi(shares=shares_SIM,log.price=log.price_SIM,
             var.soc=var.soc_SIM,log.exp=log.exp_SIM)

## changes in prices and total expenditure:
# log.price_SIM are multiplied by 1.4 between reference and current period
# log.exp_SIM are multiplied by 1.05 between reference and current period

log.price_SIM2 <- log(exp(log.price_SIM)*1.4)
log.exp_SIM2 <- log(exp(log.exp_SIM)*1.05)

## Calculation of equivalent income
equiv <- equiv.income(est,log.exp_ref=log.exp_SIM,log.exp_cur=log.exp_SIM2,
                log.price_ref=log.price_SIM,log.price_cur=log.price_SIM2)
}
}
\keyword{EASI}
\keyword{equiv}
