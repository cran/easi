residuals.easi <-
function(object=object,...){
tp <- object$summary$residuals
colnames(tp) <- object$labels.share[1:object$neq]
tp
}
