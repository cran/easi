\name{hixdata}
\alias{hixdata}
\docType{data}
\title{
Hixdata 
}
\description{
Data consist of 4,847 observations of rental-tenure single-member canadian households that had positive expenditures on rent, recreation, and transportation (For details see Arthur Lewbel and Krishna Pendakur "Tricks with Hicks : The EASI Demand System"
The American Economic Review 99(3), 827-863, 2009.
}
\usage{data(hixdata)}
\format{
  A data frame with 4847 observations on the following 26 variables.
  \describe{
    \item{\code{obs}}{number of observations.}
    \item{\code{sfoodh}}{the budget share of food at home.}
    \item{\code{sfoodr}}{the budget share of other foods.}
    \item{\code{srent}}{the budget share of rent.}
    \item{\code{soper}}{the budget share of household operations.}
    \item{\code{sfurn}}{the budget share of household furnishing and equipment.}
    \item{\code{scloth}}{the budget share of clothing.}
    \item{\code{stranop}}{the budget share of transportation operations.}
    \item{\code{srecr}}{the budget share of recreation.}
    \item{\code{spers}}{the budget share of personal care}
    \item{\code{pfoodh}}{the logarithm of the price of food at home.}
    \item{\code{pfoodr}}{the logarithm of the price of other foods.}
    \item{\code{prent}}{the logarithm of the price of rent.}
    \item{\code{poper}}{the logarithm of the price of household operations.}
    \item{\code{pfurn}}{the logarithm of the price of household furnishing and equipment.}
    \item{\code{pcloth}}{the logarithm of the price of clothing.}
    \item{\code{ptranop}}{the logarithm of the price of transportation operations.}
    \item{\code{precr}}{the logarithm of the price of recreations.}
    \item{\code{ppers}}{the logarithm of the price of personal care.}
    \item{\code{log_y}}{the logarithm of total expenditure.}
    \item{\code{age}}{the person's age minus 40.}
    \item{\code{hsex}}{the sex dummy equal to one for men.}
    \item{\code{carown}}{a dummy equal to one if real gasoline expenditures (at 1986 gasoline prices) are less than 50 dollars.}
    \item{\code{time}}{a time variable equal to the calendar year minus 1986.}
    \item{\code{tran}}{a social assistance dummy equal to one if government transfers are greater than 10 percent of gross income.}
    \item{\code{wgt}}{weighting variable.}
  }
}
%\details{
%}
\source{
Arthur Lewbel and Krishna Pendakur "Tricks with Hicks : The EASI Demand System"
The American Economic Review 99(3), 827-863, 2009
}
\references{
Arthur Lewbel and Krishna Pendakur "Tricks with Hicks : The EASI Demand System"
The American Economic Review 99(3), 827-863, 2009
}
\examples{
data(hixdata)
## maybe str(hixdata) ; plot(hixdata) ...
}
\keyword{datasets}
