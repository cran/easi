easi <-
function(shares=shares,log.price=log.price,var.soc=var.soc,log.exp=log.exp,y.power=FALSE,labels.share=FALSE,labels.soc=FALSE,py.inter=FALSE,zy.inter=FALSE,pz.inter=FALSE,interpz=FALSE){

## y.power = hightest power of y
## nsoc = mumber of demographics variables
## neq =  number of equations (without the last item)
if (!y.power){ny <- 3} else ny <- y.power
nsoc <- ncol(var.soc)
neq <- ncol(shares)-1

# Number of observations
n <- length(log.exp)

## Matrix of socio-demographic variables
z <- matrix(0,n,nsoc)
for (i in 1:nsoc)
  z[,i] <- var.soc[,i]
LABELS.Z <- c()
for (i in 1:nsoc)
LABELS.Z <-c(LABELS.Z,paste("z",i,sep=""))
colnames(z) <- LABELS.Z

## Labels or names of the budget shares: s1 - sneq
LABELS.W <- c()
for (i in 1:(neq+1))
LABELS.W <-c(LABELS.W,paste("s",i,sep=""))

LABELS.P <- c()
for (i in 1:(neq+1))
LABELS.P <-c(LABELS.P,paste("p",i,sep=""))

## Matrix of budget shares, "s", and the matrix of log.prices, "p"
s <- matrix(0,n,neq+1)
p <- matrix(0,n,neq+1)
for (i in 1:(neq+1)){
s[,i] <- shares[,i]
p[,i] <- log.price[,i]
}
colnames(s) <- LABELS.W
colnames(p) <- LABELS.P

## labels.share
if (length(labels.share)==1){
labels.share <- LABELS.W}

## labels.price
labels.price <- rep(0,neq+1)
for (i in 1:(neq+1)){
labels.price[i] <- paste("p",labels.share[i],sep="")}

## labels.soc
if (length(labels.soc)==1){
labels.soc <- LABELS.Z}

## interpz
interpz <- ifelse((length(interpz)>1),interpz,c(1:nsoc))

#####################################################################
## Convergence criteria for iterated method. If model includes     ##
## interactive terms,  we also define the the convergence variable ##
## (implicit utility or parameter)                                 ##
#####################################################################
interact <- ((zy.inter)||(py.inter)||(pz.inter))
if (interact){
conv_param <- 1
conv_y     <- 0
conv_crit  <- 0.000001
} else conv_crit  <- 0.000001

## Price normalisation: "np[,i]-np[,neq]"
np <- matrix(0,n,neq)
for (i in 1:neq)
 np[,i] <- p[,i]-p[,neq+1]

LABELS.np <- c()
for (i in 1:neq){
 LABELS.np <- c(LABELS.np,paste("np",i,sep=""))
 }
colnames(np) <- LABELS.np

## Backup of price matrix "np_backup" for future reference
np_backup <- matrix(0,n,neq)
for (i in 1:neq)
  np_backup[,i] <- np[,i]
tempo <- c()
for (i in 1:neq)
 tempo <- c(tempo,paste("np",i,"_backup",sep=""))
colnames(np_backup) <- tempo

## Initialization of matrix Ap to 0
Ap <- matrix(0,1,neq+1)
for (i in 1:(neq+1))
  Ap[i] <- 0

## Initialization of matrix  Bp and vector pBp to 0
if (interact){
Bp <- matrix(0,1,neq+1)
for (i in 1:(neq+1))
  Bp[i] <- 0
 pBp <- 0
}

## Initialization of matrix pAp to 0
pAp <- 0


########################################################################
## Creation of interaction variables "np[,i]*z[,j]"                   ##
## Not all variables need to be interacted.  Those we wish to         ##
## interact with prices will be specified by interpz=c().             ##
## For example, writing interpz=c(1,2) means that we wish to interact ##
## the first two demographic variables with the prices.               ##
########################################################################

if (pz.inter){
npz <- matrix(0,n,neq*length(interpz))
for (i in 1:neq){
  for (j in interpz){
   npz[,match(j,interpz)+(i-1)*length(interpz)]  <- np[,i]*z[,j]}}

## Creation of the names of the interaction variables "np[,i]*z[,j]"
tempo <- c()
for (i in 1:neq){
   for (j in interpz){
   tempo <- c(tempo,paste("np",i,"z",j,sep=""))}}

LABELS.npz <- tempo
colnames(npz) <- tempo

tempo <- c()
for (j in interpz){
   for (i in 1:neq){
   tempo <- c(tempo,paste("np",i,"z",j,sep=""))}}
LABELS.npz2 <- tempo
}  else npz <- LABELS.npz <- LABELS.npz2 <- c()


## Computation of  y_stone=x-p'w,
## and of the instrument, y_tilda=x-p'w^bar
y_stone <- log.exp; y_tilda <- log.exp;

mean_s <- matrix(0,n,neq+1)
for (i in 1:(neq+1))
  mean_s[,i] <- mean(s[,i])

for (i in 1:(neq+1)){
y_tilda <- y_tilda - mean_s[,i]*p[,i]
y_stone <- y_stone - s[,i]*p[,i]
}

## Creation of y^i et y_inst^i
y      <- y_stone
y_inst <- y_tilda

YY <- matrix(0,n,ny)
Yinst <- matrix(0,n,ny)
for (i in 1:ny){
YY[,i] <- y^i
Yinst[,i] <- y_inst^i
}

LABELS.YY <- c()
LABELS.Yinst <- c()
for (i in 1:ny){
 LABELS.YY <- c(LABELS.YY,paste("y",i,sep=""))
 LABELS.Yinst <- c(LABELS.Yinst,paste("y_inst",i,sep=""))
}
 
colnames(YY) <- LABELS.YY
colnames(Yinst) <- LABELS.Yinst

## Creation of y*z and z*y_inst only if required.
if (zy.inter){
yz <- matrix(0,n,nsoc)
yzinst <- matrix(0,n,nsoc)
for (i in 1:nsoc){
  yz[,i]  <- y*z[,i]
  yzinst[,i]  <- y_inst*z[,i]}

tempo <- c()
tempo2 <- c()
for (i in 1:nsoc){
  tempo <- c(tempo,paste("yz",i,sep=""))
  tempo2 <- c(tempo2,paste("yzinst",i,sep=""))
}

colnames(yz) <- tempo
colnames(yzinst) <- tempo2
} else yz <- yzinst <- c()

## Creation of  y*p and y_inst*p only if required
if (py.inter){
ynp <- matrix(0,n,neq)
for (i in 1:neq)
   ynp[,i]  <- y*np[,i]

tempo <- c()
for (i in 1:neq)
   tempo <- c(tempo,paste("ynp",i,sep=""))

LABELS.ynp <- tempo
colnames(ynp) <- tempo
 
ynpinst <- matrix(0,n,neq)
for (i in 1:neq)
   ynpinst[,i]  <- y_inst*np[,i]

tempo <- c()
for (i in 1:neq)
   tempo <- c(tempo,paste("ynpinst",i,sep=""))

colnames(ynpinst) <- tempo
} else ynp <- ynpinst <- c()

# List of covariates for each equation (labels.share)
form5 <- c("cste")
for (i in 1:ny)
  form5 <- c(form5,colnames(YY)[i])
form6 <- form5
for (i in 1:nsoc)
  form6 <- c(form6,colnames(z)[i])
form7 <- form6
if (zy.inter){
for (i in 1:nsoc)
  form7 <- c(form7,colnames(yz)[i])}
form8 <- form7
for (i in 1:neq)
  form8 <- c(form8,colnames(np)[i])
form9 <- form8
if (py.inter){
for (i in 1:neq)
  form9 <- c(form9,colnames(ynp)[i])}
form10 <- form9
if (pz.inter){
for (i in 1:(neq*length(interpz)))
  form10 <- c(form10,LABELS.npz2[i])
}

varlist <- form10

## Number of variables per equation       
dim_varlist <- length(varlist)

## Constraints of symmetry on price coefficients
TT <- matrix(0,2,neq*(neq-1)/2)
k <- 0
for (i in 1:(neq-1)){
  for (j in ((i+1):neq)){
 k <- k+1
    aa <- paste("eq",i,"_np",j,"-","eq",j,"_np",i,"=0",sep="")
TT[1,k] <- aa
}}
TT <- t(TT)
TT <- TT[,1]

## Constraints of symmetry on y*p coefficients only if required
if (py.inter){
TT2 <- matrix(0,2,neq*(neq-1)/2)
k <- 0
for (i in 1:(neq-1)){
  for (j in ((i+1):neq)){
 k <- k+1
    aa <- paste("eq",i,"_ynp",j,"-","eq",j,"_ynp",i,"=0",sep="")
TT2[1,k] <- aa
}}
TT2 <- t(TT2)
TT2 <- TT2[,1]
} else TT2 <- c()

## Constraints of symmetry on p*z coefficients only if required
if (pz.inter){
  TT3 <- matrix(0,2,neq*(length(interpz)*(neq-1))/2)
k <- 0

for (t in interpz){
for (i in 1:(neq-1)){
  for (j in ((i+1):neq)){
 k <- k+1
    aa <- paste("eq",i,"_np",j,"z",t,"-","eq",j,"_np",i,"z",t,"=0",sep="")
TT3[1,k] <- aa
}}}
TT3 <- t(TT3)
TT3 <- TT3[,1]} else TT3 <- c()

## full list of constraints
Rmat <- c(TT,TT2,TT3)

## Creation of system of equations
form1 <- c()
for (i in 1:ny)
  form1 <- paste(form1,"+",colnames(YY)[i])
form2 <- form1
for (i in 1:nsoc)
  form2 <- paste(form2,"+",colnames(z)[i])
form3 <- form2
if (zy.inter){
for (i in 1:nsoc)
  form3 <- paste(form3,"+",colnames(yz)[i])}
form4 <- form3
for (i in 1:neq)
  form4 <- paste(form4,"+",colnames(np)[i])
form5 <- form4
if (py.inter){
for (i in 1:neq)
  form5 <- paste(form5,"+",colnames(ynp)[i])}
form6 <- form5
if (pz.inter){
for (i in 1:(neq*length(interpz)))
  form6 <- paste(form6,"+",LABELS.npz2[i])
}


system <- list()
for (i in 1:neq){
  system <- c(system, list(formula(paste(paste("eqS",i,sep=""), "<-", paste("s",i,sep=""), "~", form6))))
}

## Creation of the list of instruments for the 3SLS estimation
tempo <- c()
for (i in 1:nsoc)
tempo <- c(tempo,paste("z",i,sep=""))
 colnames(z) <- tempo
form11 <- c()
for (i in 1:nsoc)
  form11 <- paste(form11,"+",colnames(z)[i])
form22 <- form11
for (i in 1:neq)
  form22 <- paste(form22,"+",colnames(np)[i])
form33 <- form22
if (pz.inter){
for (i in 1:(neq*length(interpz)))
  form33 <- paste(form33,"+",LABELS.npz2[i])}
form44 <- form33
for (i in 1:ny)
  form44 <- paste(form44,"+",colnames(Yinst)[i])
form55 <- form44
if (zy.inter){
for (i in 1:nsoc)
  form55 <- paste(form55,"+",colnames(yzinst)[i])}
form66 <- form55
if (py.inter){
for (i in 1:neq)
  form66 <- paste(form66,"+",colnames(ynpinst)[i])}

zlist <- form66
zlist <- paste("~",zlist)


## Creation of the internal database to estimate the complete system
new.data <- new.data_dep <- data.frame(cbind(s,YY,z,yz,np,ynp,npz,Yinst,yzinst,ynpinst,np_backup))

## Initialization of vectors and criteria before creating the instruments
y <- y_stone
y_backup <- y_stone
y_old <- y_stone
y_change <- 0
crit_test <- 1
iter <- 0
conv_y <- 1

## Creation of instruments (iterative linear 3sls)
cat("\n","*** Please wait during the creation of final instruments... *** ","\n")
while (crit_test>conv_crit) {
  
iter <- iter+1
  fit3sls <- systemfit( system, "3SLS", inst = formula(zlist), data = new.data,restrict.matrix=Rmat)

if (interact){  
if (iter>1) params_old <- params}

params=fit3sls$coefficients

pAp <- 0
pBp <- 0
y_old <- y


if (!interact){
#********* Predicted Values ********
pred <- predict(fit3sls,new.data)
shat <- matrix(0,n,neq)
for (i in 1:neq)
  shat[,i] <- pred[,i]

for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
}
#********* Predicted Values with p=0 and no interactions ********
pred <- predict(fit3sls,new.data)
shat_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_p0[,i] <- pred[,i]
}


if (interact){
#******** y^i = 1 *********
for (j in 1:ny)
  new.data[,paste("y",j,sep="")] <- 1}

if (py.inter){  
#********* y*p = p ******** 
for (j in 1:neq)
  new.data[,paste("ynp",j,sep="")] <- new.data[,paste("np",j,sep="")]}

if (zy.inter){
#********* y*z = z ********
for (j in 1:nsoc)
  new.data[,paste("yz",j,sep="")] <- new.data[,paste("z",j,sep="")]}

#********* Predicted Values with y=1 and interactions ********
pred <- predict(fit3sls,new.data)
shat_y1 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y1[,i] <- pred[,i]

#********* all p variables (p, p*z, p*y) are set to 0 ********
for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
if (py.inter) new.data[,paste("ynp",i,sep="")] <- 0
}

if (pz.inter){
for (i in 1:neq){
  for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- 0
}}}

#********* Predicted Values with y=1 and p variables are set to 0 and interactions********
pred <- predict(fit3sls,new.data)
shat_y1_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y1_p0[,i] <- pred[,i]

#******** p and p*z are restored  ******** p=p_backup & p*z=p_backup*z *****
for (i in 1:neq)
new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]

if (pz.inter){
for (i in 1:neq){
for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
}}}

#******** y variables and its interactions are set to 0 ********* y^i=y*p=y*z=0 ****
if (interact){
for (i in 1:ny)
new.data[,paste("y",i,sep="")] <- 0

if (py.inter){
for (i in 1:neq)
new.data[,paste("ynp",i,sep="")] <- 0}

if (zy.inter){
for (i in 1:nsoc)
new.data[,paste("yz",i,sep="")] <- 0}


#******** Predicted Values with y=0 and its interactions ********
pred <- predict(fit3sls,new.data)
shat_y0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y0[,i] <- pred[,i]

#*********** p variables are set to 0  *************
for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
new.data[,paste("ynp",i,sep="")] <- 0
}

for (i in 1:neq){
  for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- 0
}}

#******** Predicted Values with p=0 and its interactions ****
pred <- predict(fit3sls,new.data)
shat_y0_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y0_p0[,i] <- pred[,i]
}

#*********** p variables are restored : p=p_backup
for (i in 1:neq)
new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]

#********** Ap, pAp, Bp, pBp ******
if (!interact){shat_y0 <- shat; shat_y0_p0 <- shat_p0}
Ap <- matrix(0,n,neq)
for (i in 1:neq)
  Ap[,i] <- shat_y0[,i]-shat_y0_p0[,i]

pAp <- 0
#********** pAp=pAp+p*Ap ******
for (i in 1:neq)
  pAp <- pAp+new.data[,paste("np",i,sep="")]*Ap[,i]

if (interact){
Bp <- matrix(0,n,neq)
for (i in 1:neq)
  Bp[,i] <- shat_y1[,i]-shat_y1_p0[,i]-(shat_y0[,i]-shat_y0_p0[,i])

pBp <- 0
#********** pBp=pBp+p*Bp ******
for (i in 1:neq)
  pBp <- pBp+new.data[,paste("np",i,sep="")]*Bp[,i]
} else Bp <- pBp<- 0


pAp <- round(1000000*pAp+0.5)/1000000
pBp <- round(1000000*pBp+0.5)/1000000

#********* Update of y *****
y <- (y_stone+0.5*pAp)/(1-0.5*pBp)

#********* Update of y^i *****
for (i in 1:ny)
          new.data[,paste("y",i,sep="")] <- y^i

if (zy.inter){
#********* Update of y*z *****
for (i in 1:nsoc)
          new.data[,paste("yz",i,sep="")] <- y*new.data[,paste("z",i,sep="")]}

if (py.inter){
#********* y*p=y*p_backup ***********
for (i in 1:neq)
new.data[,paste("ynp",i,sep="")] <-y*new.data[,paste("np",i,"_backup",sep="")]}

if (pz.inter){
#********* z*p=z*p_backup ***********
for (i in 1:neq){
for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
}}}

if (interact){
#********* Update of crit_test if conv_param=1 *******
if ((iter>1)&(conv_param==1)){
params_change <- (params-params_old)
crit_test <- sum(params_change^2)
}}

#********* Update of y_change ********
y_change <- abs(y-y_old)

#********* Update of crit_test if conv_y=1 *******
if (conv_y==1) crit_test <- max(y_change)

#********* Update of crit_test if conv_y=1 (no interactions) *******
if (!interact) crit_test <- max(y_change)

#********* After the first iteration, conv_param replace conv_y *******
conv_y <- 0

#********* Interface => convergence *******
cat("iteration = ", iter, "\n")
cat("crit_test = ", crit_test, "\n")

#********* Creation of instruments  yinst yinst*p and yinst*z *****
y_inst <- (y_tilda+0.5*pAp)/(1-0.5*pBp)
Yinst <- matrix(0,n,ny)
for (i in 1:ny){
Yinst[,i] <- y_inst^i
}

if (py.inter){
#********* ypinst = y_inst*p **********
for (i in 1:neq)
new.data[,paste("ynpinst",i,sep="")] <- y_inst*new.data[,paste("np",i,sep="")]}

if (zy.inter){  
#********* yzinst = y_inst*z **********
for (i in 1:nsoc)
  new.data[,paste("yzinst",i,sep="")] <- y_inst*new.data[,paste("z",i,sep="")]
}

}

## Final estimate with previous instruments
new.data <- new.data_dep

#*** Initialization of vectors, parameters and criteria *****
if (interact){y_old <- y} else y_old <- y_stone
y_change <- 0
iter <- 0
crit_test <- 1

cat("\n","*** Creation of final instruments successfully completed... ***","\n")
cat("\n","*** Please wait during the estimation... ***","\n")
while (crit_test>conv_crit) {
iter <- iter+1
  fit3sls <- systemfit( system, "3SLS", inst = formula(zlist), data = new.data,restrict.matrix=Rmat)

if (interact){  
if (iter>1) params_old <- params}

params=fit3sls$coefficients

pAp <- 0
pBp <- 0
y_old <- y

if (!interact){
#********* Predicted Values ********
pred <- predict(fit3sls,new.data)
shat <- matrix(0,n,neq)
for (i in 1:neq)
  shat[,i] <- pred[,i]

for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
}

#********* Predicted Values with p=0 and no interactions ********
pred <- predict(fit3sls,new.data)
shat_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_p0[,i] <- pred[,i]
}

if (interact){
#******** y^i = 1 *********
for (j in 1:ny)
  new.data[,paste("y",j,sep="")] <- 1}

if (py.inter){  
#********* y*p = p ******** 
for (j in 1:neq)
  new.data[,paste("ynp",j,sep="")] <- new.data[,paste("np",j,sep="")]}

if (zy.inter){
#********* y*z = z ********
for (j in 1:nsoc)
  new.data[,paste("yz",j,sep="")] <- new.data[,paste("z",j,sep="")]}

#********* Predicted Values with y=1 and interactions ********
pred <- predict(fit3sls,new.data)
shat_y1 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y1[,i] <- pred[,i]

#********* all p variables (p, p*z, p*y) are set to 0 ********
for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
if (py.inter) new.data[,paste("ynp",i,sep="")] <- 0
}

if (pz.inter){
for (i in 1:neq){
  for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- 0
}}}

#********* Predicted Values with y=1 and p variables are set to 0 and interactions********
pred <- predict(fit3sls,new.data)
shat_y1_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y1_p0[,i] <- pred[,i]

#******** p and p*z are restored  ******** p=p_backup & p*z=p_backup*z *****
for (i in 1:neq)
new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]

if (pz.inter){
for (i in 1:neq){
for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
}}}

#******** y variables and its interactions are set to 0 ********* y^i=y*p=y*z=0 ****
if (interact){
for (i in 1:ny)
new.data[,paste("y",i,sep="")] <- 0

if (py.inter){
for (i in 1:neq)
new.data[,paste("ynp",i,sep="")] <- 0}

if (zy.inter){
for (i in 1:nsoc)
new.data[,paste("yz",i,sep="")] <- 0}

#******** Predicted Values with y=0 and its interactions ********
pred <- predict(fit3sls,new.data)
shat_y0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y0[,i] <- pred[,i]

#*********** p variables are set to 0  *************
for (i in 1:neq){
        new.data[,paste("np",i,sep="")] <- 0
new.data[,paste("ynp",i,sep="")] <- 0
}

for (i in 1:neq){
  for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- 0
}}

#******** Predicted Values with p=0 and its interactions ****
pred <- predict(fit3sls,new.data)
shat_y0_p0 <- matrix(0,n,neq)
for (i in 1:neq)
  shat_y0_p0[,i] <- pred[,i]
}

#*********** p variables are restored : p=p_backup
for (i in 1:neq)
new.data[,paste("np",i,sep="")] <-new.data[,paste("np",i,"_backup",sep="")]

#********** Ap, pAp, Bp, pBp ******
if (!interact){shat_y0 <- shat; shat_y0_p0 <- shat_p0}
Ap <- matrix(0,n,neq)
for (i in 1:neq)
  Ap[,i] <- shat_y0[,i]-shat_y0_p0[,i]

pAp <- 0
#********** pAp=pAp+p*Ap ******
for (i in 1:neq)
  pAp <- pAp+new.data[,paste("np",i,sep="")]*Ap[,i]

if (interact){
Bp <- matrix(0,n,neq)
for (i in 1:neq)
  Bp[,i] <- shat_y1[,i]-shat_y1_p0[,i]-(shat_y0[,i]-shat_y0_p0[,i])

pBp <- 0
#********** pBp=pBp+p*Bp ******
for (i in 1:neq)
  pBp <- pBp+new.data[,paste("np",i,sep="")]*Bp[,i]
} else Bp <- pBp<- 0


pAp <- round(1000000*pAp+0.5)/1000000
pBp <- round(1000000*pBp+0.5)/1000000

#********* Update of y *****
y <- (y_stone+0.5*pAp)/(1-0.5*pBp)

#********* Update of y^i *****
for (i in 1:ny)
          new.data[,paste("y",i,sep="")] <- y^i

if (zy.inter){
#********* Update of y*z *****
for (i in 1:nsoc)
          new.data[,paste("yz",i,sep="")] <- y*new.data[,paste("z",i,sep="")]}

if (py.inter){
#********* y*p=y*p_backup ***********
for (i in 1:neq)
new.data[,paste("ynp",i,sep="")] <-y*new.data[,paste("np",i,"_backup",sep="")]}

if (pz.inter){
#********* z*p=z*p_backup ***********
for (i in 1:neq){
for (j in interpz){
new.data[,paste("np",i,"z",j,sep="")] <- new.data[,paste("np",i,"_backup",sep="")]*new.data[,paste("z",j,sep="")]
}}}

if (interact){
#********* Update of crit_test if conv_param=1 *******
if ((iter>1)&(conv_param==1)){
params_change <- (params-params_old)
crit_test <- sum(params_change^2)
}}

#********* Update of y_change ********
y_change <- abs(y-y_old)

#********* Update of crit_test if conv_y=1 *******
if (conv_y==1) crit_test <- max(y_change)

#********* Update of crit_test if conv_y=1 (no interactions) *******
if (!interact) crit_test <- max(y_change)

#********* After the first iteration, conv_param replace conv_y *******
conv_y <- 0

#********* Interface => convergence *******
cat("iteration = ", iter, "\n")
cat("crit_test = ", crit_test, "\n")

}
cat("\n","*** Estimation successfully completed ***","\n")
#********* final estimate ****************
fit3sls <- systemfit( system, "3SLS", inst = formula(zlist), data = new.data,restrict.matrix=Rmat)

## Calculation of w_j
W=matrix(0,n,neq)
for (i in 1:neq) W[,i] <- predict(fit3sls)[,i]

colnames(W) <- labels.share[1:neq]


## Preparation of the display of results
VARS=c("Constante")
for (i in 1:ny)
VARS <- c(VARS,paste("y^",i,sep=""))
for (i in 1:nsoc)
VARS <- c(VARS,labels.soc[i])
if (zy.inter)
for (i in 1:nsoc)
  VARS <- c(VARS,paste("y*",labels.soc[i],sep=""))
for (i in 1:neq)
  VARS <- c(VARS,labels.price[i])
if (py.inter){
for (i in 1:neq)
  VARS <- c(VARS,paste("y*",labels.price[i],sep=""))}
if (pz.inter){
for (i in interpz){
  for (j in 1:neq){
  VARS <- c(VARS,paste(labels.soc[i],"*",labels.price[j],sep=""))
}}
}

VARS2 <- c("~")
for (i in 1:length(VARS))
VARS2 <- paste(VARS2,"+",VARS[i])

tempo <- c()
for (i in 1:nsoc)
tempo <- c(tempo,paste("z",i,sep=""))
 colnames(z) <- tempo
form11 <- c()
for (i in 1:nsoc)
  form11 <- c(form11,labels.soc[i])
form22 <- form11
for (i in 1:neq)
  form22 <- c(form22,labels.price[i])
form33 <- form22
if (pz.inter){
for (i in interpz){
  for (j in 1:neq){
  form33 <- c(form33,paste(labels.soc[i],"*",labels.price[j],sep=""))
}}}
form44 <- form33
for (i in 1:ny)
  form44 <- c(form44,paste("inst_y^",i,sep=""))
form55 <- form44
if (zy.inter){
for (i in 1:nsoc)
  form55 <- c(form55,paste("inst_y*",labels.soc[i],sep=""))}
form66 <- form55
if (py.inter){
for (i in 1:neq)
  form66 <- c(form66,paste("y*",labels.price[i],sep=""))}

VARINST <- form66

VARINST2 <- c("~")
for (i in 1:length(VARINST))
VARINST2 <- paste(VARINST2,"+",VARINST[i])

a=summary(fit3sls)
for (i in 1:neq){
a$eq[[i]][1] <- labels.share[i]
rownames(a$eq[[i]][8][[1]]) <- VARS
a$eq[[i]][3] <- VARS2
a$eq[[i]][4] <- VARINST2
}

Residuals=summary(fit3sls)$residuals
colnames(Residuals) <- labels.share[1:neq]



##########################################################################
# The EASI function returns "Result", an object of class "easi", namely  #
# a list of results that includes:                                       # 
# - Residuals                                                            #
# - The covariance matrix of estimated parameters                        #
# - The matrix of the fitted budget shares                               #
# - a summary of the result (a)				                 #
# - the implicit utility (y)					         #
# - the list of all variables                                            #
# - others objects used in the estimation : var.soc (demographic varia-  #
# bles), shares (matrix of budget shares), log.price (matrix of prices),     #
# neq (number of equations without the last equation), y.power (hightest      #
# power of y), nsoc (number of demographic variables), interact (a dum-  #
# my equal to one if there are interaction), py.inter (a dummy equal to  #
# one if there are interactions between Prices and y, pz.inter (a dummy  #
# equal to one if there are interactions between Prices and demographic  #
# variables, zy.inter (a dummy equal to one if there are interactions    #
# between demographics and y, interpz (that allows to choose the varia-  #
# bles Z to crossed with the price), fit3sls (object of class systemfit) #
# log.exp (the logaritm of total expenditure), labels.price (the names of      #
# prices), labels.soc (the names of demographic variables), labels (the  #
# names of budget shares), dim_varlist (nombre de variables)             #
##########################################################################

Result <-  list(
Residuals=Residuals,
CoefCov=summary(fit3sls)$coefCov,
fitted.w=W,
summary=a,
y=y,
varlist=varlist,
var.soc=var.soc,
shares=shares,
log.price=log.price,
neq=neq,
y.power=ny,
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
VARS=VARS
)
class(Result) <- c("easi","systemfit")
Result

}
