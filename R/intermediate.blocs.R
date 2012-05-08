intermediate.blocs <-
function(object=object,log.price=log.price,var.soc=var.soc,log.exp=log.exp){

shares <- object$shares
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
log.exp <- object$log.exp
labels.price <- object$labels.price
labels.soc <- object$labels.soc
labels.share <- object$labels.share
dim_varlist <- object$dim_varlist
y <- object$y

n <- length(log.exp)

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


## Recovery of coefficients for the variables p and p * z Note: the first element of Z is a constant to capture the direct price effects

my.array <-array (0, dim = c (nsoc+1,(neq+1),(neq+1)))

## Step 1: Recovery of coefficients of p variables
a0 <- matrix(0,neq,neq)
for (i in 1:neq){
  for (j in 1:neq){
a0[i,j] <- coef[paste("eq",i,"_np",j,sep="")]
my.array[,,i][1,j] <- a0[i,j]
}}
my.array[,,neq+1][1,1:neq] <- 0-apply(a0,2,sum)

for (i in 1:neq){
my.array[,,i][1,neq+1] <- my.array[,,neq+1][1,i] 
}
my.array[,,neq+1][1,neq+1] <- 0-sum(my.array[,,neq+1][1,1:neq])



## Step 2: Recovery of coefficients of p*z variables only if required
if (pz.inter){
for (i in 1:neq){
  for (j in interpz){
    for (k in 1:neq){
    my.array[,,i][j+1,k] <- coef[paste("eq",i,"_np",k,"z",j,sep="")]
  }}}
}

for (t in interpz){
for (i in 1:neq){
for (j in 1:neq){
my.array[,,neq+1][t+1,i] <- my.array[,,neq+1][t+1,i]-my.array[,,j][t+1,i]
}
}
}

for (t in interpz){
for (i in 1:neq){
my.array[,,i][t+1,neq+1] <- my.array[,,neq+1][t+1,i]
}
}

for (t in interpz){
for (i in 1:neq){
my.array[,,neq+1][t+1,neq+1] <- my.array[,,neq+1][t+1,neq+1]-my.array[,,i][t+1,neq+1]
}
}

## construction of the sum "sum_j sum_k sum_t a_jkt z_t p_j p_k" (calculation of y)
## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
a <- my.array
tot = 0
for (j in 1:neq){
for (k in 1:neq){
for (t in c(1:(nsoc+1))){
tempo <- a[t,k,j]*P[,k]*P[,j]*Z[,t]
tot <- tot+tempo
}}}

## Recovery of coefficients of p*y variables only if required
bjk=matrix(0,neq+1,neq+1)
tot2 = 0
if (py.inter){
for (i in 1:neq){
  for (j in 1:neq){
   bjk[j,i] <- coef[paste("eq",i,"_ynp",j,sep="")]
 }}

  for (j in 1:neq){
   bjk[j,neq+1] <- 0-sum(bjk[j,1:neq])
 }
 
  for (j in 1:neq){
   bjk[neq+1,j] <- bjk[j,neq+1]
 }

   bjk[neq+1,neq+1] <- 0-sum(bjk[neq+1,1:neq])
 
colnames(bjk)<- c(noms,"Others")

## construction of the sum "sum_j sum_k  b_jk p_j p_k" (calculation of y)
## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
for (j in 1:neq){
for (k in 1:neq){
tempo <- bjk[j,k]*P[,j]*P[,k]
tot2 <- tot2+tempo
}}
}

## construction of the sum "sum_j w_j p_j" (calculation of y)
## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
tot0 = 0
for (j in 1:neq){
tempo <- w[,j]*P[,j]
tot0 <- tot0+tempo
}

## Calculation of y
## "EASI made EASIER" (PENDAKUR 2008 -  page 11 formula 22)
if (interact){
y <- (lnx - tot0+1/2*tot)/(1-1/2*tot2)} else y <- (lnx-tot0+1/2*tot)

## Recovery of coefficients of y^r variables (calculation of w_j)
bjr=matrix(0,y.power,neq+1)
for (i in 1:neq){
  for (j in 1:y.power){
   bjr[j,i] <- coef[paste("eq",i,"_y",j,sep="")]
 }}

  for (j in 1:y.power){
   bjr[j,neq+1] <- 0-sum(bjr[j,1:neq])
 }


colnames(bjr)<- c(noms,"others")

## Recovery of coefficients of z variables (calculation of w_j)
gjt=matrix(0,nsoc,neq+1)
for (i in 1:neq){
  for (j in 1:nsoc){
   gjt[j,i] <- coef[paste("eq",i,"_z",j,sep="")]
 }}

  for (j in 1:nsoc){
   gjt[j,neq+1] <- 0-sum(gjt[j,1:neq])
 }

colnames(gjt)<- c(noms,"others")


## Recovery of coefficients of y*z variables (calculation of w_j) only if required
hjt=matrix(0,nsoc,neq+1)
if (zy.inter){
for (i in 1:neq){
  for (j in 1:nsoc){
   hjt[j,i] <- coef[paste("eq",i,"_yz",j,sep="")]
 }}

 for (j in 1:nsoc){
   hjt[j,neq+1] <- 0-sum(hjt[j,1:neq])
 }

}
colnames(hjt) <- c(noms,"others")

## Recovery of the constants (calculation of w_j)
cc=c()
for (i in 1:neq)
cc=cbind(cc,coef[paste("eq",i,"_(Intercept)",sep="")])

cc=cbind(cc,1-sum(cc))

colnames(cc) <- c(noms,"others")


Result <-  list(
CoefCov=fit3sls$coefCov,
a=a,
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
fit3sls=fit3sls,
log.exp=log.exp,
labels.price=labels.price,
labels.soc=labels.soc,
labels.share=labels.share,
dim_varlist=dim_varlist,
n=n,
coef=coef,
my.array=my.array,
tot=tot,
tot2=tot2,
tot0=tot0,
bjk=bjk,
P=P,
w=w,
Z=Z,
bjr=bjr,
gjt=gjt,
hjt=hjt,
cc=cc,
lnx=lnx

)

Result

}
