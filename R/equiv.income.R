equiv.income <-
function(object=object,log.exp_ref=log.exp_ref,log.exp_cur=log.exp_cur,log.price_ref=log.price_ref,log.price_cur=log.price_cur){

## recovery of results of previous estimation
fit3sls <- object$fit3sls
varlist <- object$varlist
neq <- object$neq
y.power <- object$y.power
nsoc <- object$nsoc
interact <- object$interact
py.inter <- object$py.inter
zy.inter <- object$zy.inter
pz.inter <- object$pz.inter
interpz <- object$interpz


## Definition of ...
var.soc <- object$var.soc
shares <- object$shares
n <- length(log.exp_cur)

temp <- intermediate.blocs(object,log.price=log.price_cur,var.soc=var.soc,log.exp=log.exp_cur)
tot <- temp$tot
tot2 <- temp$tot2
tot0 <- temp$tot0
if (interact){
y <- (log.exp_cur - tot0+1/2*tot)/(1-1/2*tot2)} else y <- (log.exp_cur -tot0+1/2*tot)
A <- y

################################
################################
#####   Reference Situation  ###
################################
################################

temp <- intermediate.blocs(object,log.price=log.price_ref,var.soc=var.soc,log.exp=log.exp_ref)
tot <- temp$tot
tot2 <- temp$tot2
tot0 <- temp$tot0

#######################################
#######################################
##  calculation of equivalent income ##
#######################################
#######################################

B <- tot0-1/2*tot
C <- 1-1/2*tot2
equiv_income <- exp(A*C+B)

cat("\n","Info_1: The average Equivalent income is equal to :",mean(equiv_income),"\n")

## Validation 0: implicit utility with reference income in the reference situation (start situation)
if (interact){
y <- (log.exp_ref - tot0+1/2*tot)/(1-1/2*tot2)} else y <- (log.exp_ref -tot0+1/2*tot)

cat("\n","Info_2: The average implicit utility with reference income in","\n") 
cat("\n","           the reference situation is equal to :",mean(y),"\n")
cat("\n","*** it is the implicit utility before any changes ***","\n")

## Current Implicit utility 
cat("\n","Info_3: The Current Implicit Utility is equal to:",mean(A),"\n")

## Validation 1: implicit utility with income equivalent in the reference situation = implicit utility with contemporary income in the contemporary situation
if (interact){
y <- (log(equiv_income) - tot0+1/2*tot)/(1-1/2*tot2)} else y <- (log(equiv_income) -tot0+1/2*tot)

cat("\n","Info_4: The average implicit utility with income equivalent in","\n")
cat("\n","           the reference situation is equal to :",mean(y),"\n")
cat("\n","*** it should be equal to implicit utility in the current situation above ***","\n")


## Validation 2: implicit utility with contemporary income in the reference situation < or > implicit utility with contemporary income in the contemporary situation

if (interact){
y2 <- (log.exp_cur - tot0+1/2*tot)/(1-1/2*tot2)} else y2 <- (log.exp_cur -tot0+1/2*tot)


cat("\n","Info_5: The average implicit utility with contemporary income in","\n") 
cat("\n","           the reference situation is equal to :",mean(y2),"\n")

equiv_income
}
