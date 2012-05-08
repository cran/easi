\name{engel}
\alias{engel}
\title{
Engel curves of EASI demand system
}
\description{
The engel function uses estimates of the easi function in order to calculate and draw the Engel curves.
}
\usage{
engel(object = object, file = file, sd = FALSE, limY = FALSE)
}
\arguments{
  \item{object}{
an object of type easi which contains all the information about the previous estimation. Here, in particular, the estimates of the easi system.
}
  \item{file}{
filename with which Engel Curves ("file.pdf") are to be stored. For example, if one specifies "file=Myfile1", in the parent directory, the file "Myfile1.pdf" - which contains Engel curves - is generated.
}
  \item{sd}{
sd (for standard deviations) must be specified as TRUE if we want calculate the standard deviations (FALSE otherwise).
}
  \item{limY}{
the  limits of the Y-axis of the graphical representations of Engel curves. For example, limY=c(0,0.5,0.2,0.5,0,1) means that the first Engel curve will be represented on a graph whose limits on the Y-axis are [0;0.5], the second Engel curve will be represented on a graph whose limits on the Y-axis are [0.2;0.5], and so on. By default, all Y-axis limits are [0;third quartile of observed budget shares]. It is also a good way to evaluate the goodness of fit.
}
}
\value{
The matrix of fitted budget shares.
}
\details{ 
The equations of the Engel curves are presented in Pendakur and LewbeL. 
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


#******* engel curves *************************
eng1 <- engel(est,file="graph_engels_curves",sd=FALSE)
}
}
\keyword{EASI}
\keyword{Engel}
