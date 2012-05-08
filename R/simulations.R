simulations <-
function(object=object,log.price_new=log.price_new,var.soc_new=var.soc_new,log.exp_new=log.exp_new){

fit3sls <- object$fit3sls
varlist <- object$varlist
shares <- object$shares
neq <- object$neq
y.power <- object$y.power
nsoc <- object$nsoc
interact <- object$interact
py.inter <- object$py.inter
zy.inter <- object$zy.inter
pz.inter <- object$pz.inter
interpz <- object$interpz
labels.price <- object$labels.price
labels.soc <- object$labels.soc
labels.share <- object$labels.share
dim_varlist <- object$dim_varlist
log.exp <- object$log.exp

n <- length(log.exp)

log.price=log.price_new 
var.soc=var.soc_new
log.exp <- log.exp_new 



temp <- intermediate.blocs(object,log.price=log.price,var.soc=var.soc,log.exp=log.exp)
my.array <- temp$my.array
tot <- temp$tot
tot2 <- temp$tot2
tot0 <- temp$tot0
bjk <- temp$bjk
P <- temp$P
w <- temp$w
Z <- temp$Z
bjr <- temp$bjr
gjt <- temp$gjt
hjt <- temp$hjt
cc <- temp$cc
noms <- object$noms
lnx <- object$log.exp
y <- object$log.exp



## Calculation of w_j
W=matrix(0,n,neq)

ajk <- my.array
for (i in 1:neq){

  tot3 <- tot4 <- tot5 <- tot6 <- tot7 <- 0
  
  for (j in 1:y.power){
    tempo <- bjr[j,i]*y^j
    tot3 <- tot3+tempo}

  for (j in 1:nsoc){
    tempo <- gjt[j,i]*Z[,j+1]
    tot4 <- tot4+tempo}

  if (zy.inter){
   for (j in 1:nsoc){
    tempo <- hjt[j,i]*Z[,j+1]*y
    tot5 <- tot5+tempo}}
  
if (pz.inter){
  for (k in 1:neq){
    for (t in (1:(nsoc+1))){
      tempo <- ajk[t,k,i]*Z[,t]*P[,k]
      tot6 <- tot6+tempo
    }}} 

  if (py.inter){
  for (k in 1:neq){
    tempo <- bjk[k,i]*P[,k]*y
    tot7 <- tot7+tempo }}

 
  W[,i] <- cc[i]+tot3+tot4+tot5+tot6+tot7

}


colnames(W) <- labels.share[1:neq]

Result <-  list(
CoefCov=fit3sls$coefCov,
fit3sls=fit3sls,
y=y,
varlist=varlist,
var.soc=var.soc,
shares=shares,
log.price=log.price,
neq=neq,
y.power=y.power,
nsoc=nsoc,
interact=interact,
py.inter=py.inter,
zy.inter=zy.inter,
pz.inter=pz.inter,
interpz=interpz,
log.exp=log.exp,
labels.price=labels.price,
labels.soc=labels.soc,
labels.share=labels.share,
dim_varlist=dim_varlist,
fitted.w=W
)

Result

}
