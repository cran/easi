coef.easi <-
function(object=object,...){
neq <- object$neq
VARS <- object$VARS
VARS2 <- c()
for (i in 1:neq)
VARS2 <- c(VARS2,paste(paste("eq",i,sep=""),VARS,sep="_"))
tp <- object$summary$coefficients
rownames(tp) <- VARS2
tp
}
