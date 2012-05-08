engel <-
function(object=object,file=file,sd=FALSE,limY=FALSE){

WDELTA <- ifelse(sd,TRUE,FALSE) 
fit3sls <- object$fit3sls
varlist <- object$varlist
var.soc <- object$var.soc
shares <- object$shares
log.price <- object$log.price
neq <- object$neq
y.power <- object$y.power
nsoc <- object$nsoc
interact <- object$interact
py.inter <- object$py.inter
zy.inter <- object$zy.inter
pz.inter <- object$pz.inter
interpz <- object$interpz
log.exp <- object$log.exp
labels.price <- object$labels.price
labels.soc <- object$labels.soc
labels.share <- object$labels.share
dim_varlist <- object$dim_varlist
y <- object$y

n <- length(log.exp)


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
noms <- temp$noms
lnx <- temp$log.exp
y <- temp$y

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


## Calculation of standard deviations of the fitted budget shares (if WDELTA=TRUE) - (delta method)
if (WDELTA){
nb <- dim_varlist
MAT <- rep(1,n)
for (i in 1:y.power)
MAT <- cbind(MAT,y^i)
for (i in 1:nsoc)
MAT <- cbind(MAT,Z[,i+1])
if (zy.inter){
for (i in 1:nsoc)
MAT <- cbind(MAT,y*Z[,i+1])
}
for (i in 1:neq)
MAT <- cbind(MAT,P[,i])
if (py.inter){
for (i in 1:neq)
MAT <- cbind(MAT,y*P[,i])
}
if (pz.inter){
for (i in interpz){
  for (j in 1:neq){
MAT <- cbind(MAT,Z[,i+1]*P[,j])}}
}

nn <- 1
W_ecart <- matrix(0,n,neq)

for (i in 1:neq){
DD <- summary(fit3sls)$coefCov[nn:(nn+nb-1),nn:(nn+nb-1)]

W_e <- MAT%*%DD%*%t(MAT)

W_ecart[,i] <- sqrt(diag(W_e))
rm(W_e)

nn <- nn+nb
}
}


## Labels of W matrix
colnames(W)=noms

## Engel Curves
quant=quantile(lnx,seq(0,1,0.01))

ee=rep(1,n)
for (j in 1:n) {
for (i in 1:100)
if (lnx[j]>quant[i]) 
ee[j]=i
}

Wm=matrix(0,100,neq)
for (i in 1:100){
for (j in 1:neq)
Wm[i,j]=median(W[ee==i,j])
}

Wm_autres=1-apply(Wm,1,sum)
Wm=cbind(Wm,Wm_autres)


## Calculation of confidence intervals for fitted budget shares (if WDELTA=TRUE)
if (WDELTA){
Wme=matrix(0,100,neq+1)
for (i in 1:100){
for (j in 1:neq)
Wme[i,j]=median(W_ecart[ee==i,j])
}

for (i in 1:100){
Wme[i,neq+1]=sqrt(sum(Wme[i,1:neq]^2))
}

Wmep <- Wm+1.96*Wme
Wmem <- Wm-1.96*Wme

Wmep <- Wmep[(1:20)*5,]
Wmem <- Wmem[(1:20)*5,]
}

### management of labels.share
if (length(labels.share)<2) labels.share <- noms

limYY <- c()
if (length(limY)<2){
    for (i in 1:neq) 
      limYY <- c(limYY,c(0,summary(w[,i])[5]))} else {
      limYY <- limY}

ss <- seq(1,neq*2,by=2)

## Export of Engel curves in the parent folder under the name "file". pdf
## File name is entered on the command line

pdf(paste("./",file,".pdf"))
for(i in 1:neq)
{
# smoothing cubic
 sp <- smooth.spline(c(1:100),Wm[,i], spar = 0.9)
y.loess <- loess(Wm[,i] ~ c(1:100), span=0.75, data.frame(xxx=c(1:100), yyy=Wm[,i]))
y.predict <- predict(y.loess, data.frame(xxx=c(1:100)))
 plot(c(1:100),Wm[,i],xlab="Percentiles of total expenditure",ylab="Budget shares",col="green",ylim=c(limYY[ss[i]],limYY[ss[i]+1]))
title(main=labels.share[i])
xx <- seq(1,100, len=20)
# plot of the adjustment curve 
lines(predict(sp,xx), col = "red")
 lines(c(1:100),y.predict,col="blue")
 lines(ksmooth(c(1:100),Wm[,i], "normal", bandwidth=10), col="black")
if (WDELTA){  points(c((1:20)*5),Wmep[,i],pch="+",cex=1,col="violet")
  points(c((1:20)*5),Wmem[,i],pch="+",cex=1,col="violet")}
       }
sp <- smooth.spline(c(1:100),Wm[,i+1], spar = 0.9)
plot(c(1:100),Wm[,i+1],xlab="Percentiles of total expenditure",ylab="Budget shares",col="green")
title(main ="Others Goods")
xx <- seq(1,100, len=20)
y.loess <- loess(Wm[,i+1] ~ c(1:100), span=0.75, data.frame(xxx=c(1:100), yyy=Wm[,i]))
y.predict <- predict(y.loess, data.frame(xxx=c(1:100)))
# plot of the adjustment curve
lines(predict(sp,xx), col = "red")
lines(c(1:100),y.predict,col="blue")
lines(ksmooth(c(1:100),Wm[,i+1], "normal", bandwidth=10), col="black")
if (WDELTA){  points(c((1:20)*5),Wmep[,i+1],pch="+",cex=1,col="violet")
  points(c((1:20)*5),Wmem[,i+1],pch="+",cex=1,col="violet")}
dev.off()


## refreshment
if (WDELTA){rm(MAT); rm(W_ecart); rm(DD)}

result <- W
result
}
