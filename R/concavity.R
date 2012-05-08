concavity <-
function(object=object){

cat("\n","*** Please wait during the test of local concavity *** ","\n")

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
noms <- object$noms
lnx <- object$log.exp
y <- object$y

#** Recovery of the estimated coefficients***
Estimates <- summary(fit3sls)$coefficients

coef=Estimates[,1]
n <- length(log.exp)

## Labels or names of the equations:
noms <- c()
for (i in 1:neq)
  noms <- c(noms, paste("eq",i,sep=""))

#***** new price Matrix *******
P=log.price

#****  new sociodemographic matrix *****
Z=cbind(1,var.soc)

#**** logarithm of total expenditure *******
lnx <- log.exp

#**** Budget shares matrix *******
w=matrix(0,n,neq+1)
for (i in 1:(neq+1))
  w[,i] <- shares[,i]

## Calculation of first and second derivative of y with respect to p_j and p_k
## uprime
h1 <- my.array
toth1 <- matrix(0,n,neq)
u1p <- matrix(0,n,neq)
for (j in 1:neq){
toth1[,j] = 0
for (k in 1:neq){
for (t in c(1:(nsoc+1))){
tempo <- h1[t,k,j]*P[,k]/exp(P[,j])*Z[,t]
toth1[,j] <- toth1[,j]+tempo
}
}
u1p[,j]=-w[,j]/exp(P[,j])+toth1[,j]
}

## usecond
toth2 <- array(0,c(neq,neq,n))
u2p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
toth2[j,k,] = 0
for (t in c(1:(nsoc+1))){
tempo <- h1[t,k,j]/(exp(P[,k])*exp(P[,j]))*Z[,t]
toth2[j,k,] <- toth2[j,k,]+tempo
}
}
}
u2p=toth2

## vprime
toth3 <- matrix(0,n,neq)
v1p <- matrix(0,n,neq)
for (k in 1:neq){
tempo <- bjk[j,k]/exp(P[,j])*P[,k]
toth3[,k] <- toth3[,k]+tempo
}
for (j in 1:neq){
v1p[,j]=-toth3[,j]}


## vsecond
toth4 <- array(0,c(neq,neq,n))
v2p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
toth4[j,k,] = bjk[j,k]/(exp(P[,j])*exp(P[,k]))
}
}
v2p=toth4

S=y+tot0-1/2*tot-1/2*tot2*y

u=lnx - tot0+1/2*tot
v=1-1/2*tot2

w <- array(0,c(neq,neq,n))
wp <- array(0,c(neq,neq,n))
lp <- array(0,c(neq,neq,n))
l <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
w[j,k,] <- u1p[,j]*v-v1p[,j]*u
l[j,k,] <- v*v
wp[j,k,] <- u2p[j,k,]*v + u1p[,j]*v1p[,j]-v1p[,j]*u1p[,j]-v2p[j,k,]*u  
lp[j,k,] <- 2*v1p[,j]*v
}
}

Hy1 <- wp/l
Hy2 <- (wp*l-lp*w)/(l)^2

## Calculation of of first and second derivative of S0 with respect to p_j and p_k :
S0 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S0[j,k,] <- tot0
}
}

S01p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S01p[j,k,] <- shares[,j]/exp(P[,j])
}
}

S02p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S02p[j,k,] <- -(j==k)*shares[,j]/exp(P[,j])^2
}
}

## Calculation of of first and second derivative of S1 with respect to p_j and p_k :
S1 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S1[j,k,] <- -1/2*tot
}
}

S11p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S11p[j,k,] <- -toth1[,j]
}
}

S12p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S12p[j,k,] <- -toth2[j,k,]
}
}

## Calculation of of first and second derivative of S2 with respect to p_j and p_k :
S2 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S2[j,k,] <- -1/2*tot2
}
}

S21p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S21p[j,k,] <- -toth3[,j]
}
}

S22p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S22p[j,k,] <- -toth4[j,k,]
}
}

## First derivate of S in respect with p_j
SSp <- Hy1 + S01p + S11p + S21p*y + S2*Hy1

## Second derivate of S
SS2p <- Hy2 + S02p + S12p + S22p*y + S21p*Hy1 + S21p*Hy1+S2*Hy2

w <- shares

## Calculation of first and second derivative of y with respect to p_j and p_k
## uprime
h1 <- my.array
toth1 <- matrix(0,n,neq)
u1p <- matrix(0,n,neq)
for (j in 1:neq){
toth1[,k] = 0
for (k in 1:neq){
for (t in c(1:(nsoc+1))){
tempo <- h1[t,k,j]*P[,k]/exp(P[,k])*Z[,t]
toth1[,k] <- toth1[,k]+tempo
}
}
u1p[,k]=-w[,k]/exp(P[,k])+toth1[,k]
}

## usecond
toth2 <- array(0,c(neq,neq,n))
u2p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
toth2[,k,j] = 0
for (t in c(1:(nsoc+1))){
tempo <- h1[t,k,j]/(exp(P[,j])*exp(P[,k]))*Z[,t]
toth2[k,j,] <- toth2[k,j,]+tempo
}
}
}
u2p=toth2

## vprime
toth3 <- matrix(0,n,neq)
v1p <- matrix(0,n,neq)
for (j in 1:neq){
tempo <- bjk[j,k]/exp(P[,k])*P[,j]
toth3[,j] <- toth3[,j]+tempo
}
for (k in 1:neq){
v1p[,k]=-toth3[,k]}

## vsecond
toth4 <- array(0,c(neq,neq,n))
v2p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
toth4[k,j,] = bjk[j,k]/(exp(P[,k])*exp(P[,j]))
}
}
v2p=toth4

S=y+tot0-1/2*tot-1/2*tot2*y

u=lnx - tot0+1/2*tot
v=1-1/2*tot2

wp <- array(0,c(neq,neq,n))
l <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
l[k,j,] <- v*v
wp[k,j,] <- u2p[k,j,]*v + u1p[,k]*v1p[,k]-v1p[,k]*u1p[,k]-v2p[k,j,]*u  

}
}

Hy1 <- wp/l

## Calculation of of first and second derivative of S0 with respect to p_j and p_k :
S0 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S0[k,j,] <- tot0
}
}

S01p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S01p[k,j,] <- shares[,k]/exp(P[,k])
}
}

## Calculation of of first and second derivative of S1 with respect to p_j and p_k :
S1 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S1[k,j,] <- -1/2*tot
}
}

S11p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S11p[k,j,] <- -toth1[,k]
}
}

## Calculation of of first and second derivative of S2 with respect to p_j and p_k :
S2 <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S2[k,j,] <- -1/2*tot2
}
}

S21p <- array(0,c(neq,neq,n))
for (j in 1:neq){
for (k in 1:neq){
S21p[k,j,] <- -toth3[,k]
}
}

## First derivate of S in respect with p_k
SSpp <- Hy1 + S01p + S11p + S21p*y + S2*Hy1



hess <- exp(S)*(SS2p+SSp*SSpp)

testconcav <- matrix(0,n,1)
for (i in 1:n){
testconcav[i] <- semidefiniteness(hess[,,i],positive=TRUE)
if (!testconcav[i]) testconcav[i] <- 1 else testconcav[i] <- 0
}


for (i in (10:1)/10){ 
if (mean(testconcav)<=i)
a <- paste(" Result: The cost function is concave on more than ",(i-0.1)*100,"% of the sample. ",sep="")
}

cat("\n")
a
}





