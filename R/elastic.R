elastic <-
function(object=object,type = c("price", "income", "demographics"),sd=FALSE){


  type <- match.arg(type)

EPDELTA <- ifelse(((type=="price")&(sd==TRUE)),TRUE,FALSE)
ERDELTA <- ifelse(((type=="income")&(sd==TRUE)),TRUE,FALSE)
EZDELTA <- ifelse(((type=="demographics")&(sd==TRUE)),TRUE,FALSE)


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


if (type == "price"){
### Calculation of log.price elasticities

#*** semi-elasticities with respect to log.prices
## page 13 formula 23 EASI made EASIER (Pendakur 2008)
EP=matrix(0,neq+1,neq+1)
a <- my.array 
for(i in 1:(neq+1)){
  for (k in 1:(neq+1)){
    tot10 = 0
    for (t in (1:(nsoc+1))){
      tempo <- a[t,k,i]*Z[,t]
      tot10 <- tot10+tempo
    }
      tot10 <- tot10+bjk[k,i]*y
  EP[k,i]=mean(tot10)
}}


colnames(EP) <- rownames(EP) <- c(labels.price[1:neq],"pothers")


## Calculation of standard deviations of log.price elasticities
## (Delta method) if EPDELTA=TRUE 
if (EPDELTA){
ttt <- colnames(summary(fit3sls)$coefCov)
EP_SE=matrix(0,neq+1,neq+1)
ELASTPRICE_SE=matrix(0,neq+1,neq+1)
for (i in 1:neq){
  for (j in 1:neq){
 tt <- paste("eq",i,"_np",j,sep="")

if (pz.inter){
 for (t in interpz){
  tt <- c(tt,paste("eq",i,"_np",j,"z",t,sep=""))
}}

if (py.inter){
 tt <- c(tt,paste("eq",i,"_ynp",j,sep=""))
}


 tnum <- match(tt,ttt)

  DD <- summary(fit3sls)$coefCov[tnum,tnum]

 MAT <- Z[,1]
 if (pz.inter)
 MAT <- cbind(MAT,Z[,interpz+1])
 if (py.inter)
 MAT <- cbind(MAT,y)

 
   EP_SE[i,j] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))
ELASTPRICE_SE[i,j] <- EP_SE[i,j]/mean(shares[,j])
rm(DD)}}


EP_SE[neq+1,1:neq] <- sqrt(apply(EP_SE[1:neq,1:neq]^2,1,sum))
EP_SE[,neq+1] <- EP_SE[neq+1,]
EP_SE[neq+1,neq+1] <- sqrt(sum(EP_SE[neq+1,1:neq]))

for (i in (1:(neq+1))){
ELASTPRICE_SE[i,neq+1] <- EP_SE[i,neq+1]/mean(shares[,neq+1])}

ELASTPRICE_SE[neq+1,1:neq] <- ELASTPRICE_SE[1:neq,neq+1]
ELASTPRICE_SE[neq+1,neq+1] <- sqrt(sum(ELASTPRICE_SE[1:neq,neq+1]))


colnames(ELASTPRICE_SE) <- rownames(ELASTPRICE_SE) <- colnames(EP_SE) <- rownames(EP_SE) <- c(labels.price[1:neq],"pothers")

}



#*** Normalised Slutsky matrix
#** matrix of compensated quantity derivatives with respect to unlogged log.prices)
#** own-log.price Slutsky terms in Pendakur
## page 849  "Tricks with Hicks : The EASI demand system" (Lewbel & Pendakur 2008)
EPS=EP+apply(shares[,1:(neq+1)],2,mean)%*%t(apply(shares[,1:(neq+1)],2,mean))-matrix(diag(apply(shares[,1:(neq+1)],2,mean)),neq+1,neq+1)
colnames(EPS) <- rownames(EPS) <- c(labels.price[1:neq],"pothers")


#*** Compensated (good-specific) expenditures elasticities with respect to log.prices
#** own-log.price Quant elast in Pendakur
## page 849 "Tricks with Hicks : The EASI demand system" (Lewbel & Pendakur 2008)
EPQ=solve(diag(apply(shares[,1:(neq+1)],2,mean)))%*%(EP+apply(shares[,1:(neq+1)],2,mean)%*%t(apply(shares[,1:(neq+1)],2,mean)))
colnames(EPQ) <- rownames(EPQ) <- c(labels.price[1:neq],"pothers")


#*** calculation of elasticity of good j with respect to the log.price of good i
#*** calculation of elasticity of good j with respect to income
ajk <- my.array
ELASTPRICE <- matrix(0,neq+1,neq+1)
for (i in 1:(neq+1)){
for (q in 1:(neq+1)){

C <- 0
  for (j in 1:y.power){
    tempo <- j*bjr[j,i]*y^{j-1}
    C <- C+tempo}

D <- 0
  if (zy.inter){
   for (j in 1:nsoc){
    tempo <- hjt[j,i]*Z[,j+1]
    D <- D+tempo}}

E <- 0
    for (t in (1:(nsoc+1))){
      tempo <- ajk[t,q,i]*Z[,t]
      E <- E+tempo
    } 

G <- 0
  if (py.inter){
  for (k in 1:(neq+1)){
    tempo <- bjk[k,i]*P[,k]
    G <- G+tempo }

F <- bjk[q,i]*y
}

U <- 0
 for (d in 1:(neq+1)){
    for (t in (1:(nsoc+1))){
      tempo <- ajk[t,d,i]*Z[,t]*P[,d]
      U <- U+tempo
    } }

B <- -mean(shares[,q]+U)/mean(1-1/2*tot2)-mean(y)/mean(1-1/2*tot2)*mean(G)

H <- B*(C+D+G)+E+F

ELASTPRICE[q,i] <- mean(H) / mean(shares[,i]) - as.numeric((i==q))

}
}

colnames(ELASTPRICE) <- rownames(ELASTPRICE) <- c(labels.price[1:neq],"pothers")

}










if (type=="income"){

#*** calculation of elasticity of good j with respect to the log.price of good i
#*** calculation of elasticity of good j with respect to income
ajk <- my.array
ELASTINCOME <- matrix(0,1,neq+1)
for (i in 1:(neq+1)){
for (q in 1:(neq+1)){

C <- 0
  for (j in 1:y.power){
    tempo <- j*bjr[j,i]*y^{j-1}
    C <- C+tempo}

D <- 0
  if (zy.inter){
   for (j in 1:nsoc){
    tempo <- hjt[j,i]*Z[,j+1]
    D <- D+tempo}}

E <- 0
    for (t in (1:(nsoc+1))){
      tempo <- ajk[t,q,i]*Z[,t]
      E <- E+tempo
    } 

G <- 0
  if (py.inter){
  for (k in 1:(neq+1)){
    tempo <- bjk[k,i]*P[,k]
    G <- G+tempo }

F <- bjk[q,i]*y
}

U <- 0
 for (d in 1:(neq+1)){
    for (t in (1:(nsoc+1))){
      tempo <- ajk[t,d,i]*Z[,t]*P[,d]
      U <- U+tempo
    } }

B <- -mean(shares[,q]+U)/mean(1-1/2*tot2)-mean(y)/mean(1-1/2*tot2)*mean(G)

H <- B*(C+D+G)+E+F

ELASTINCOME[1,i] <- 1 + (mean(C+D+G))/mean(shares[,i])

}
}

colnames(ELASTINCOME) <- c(labels.share[1:neq],"others")



## Calculation of income elasticities of budget shares
## page 13 formula 23 "EASI made EASIER" (Pendakur 2008)
ER=matrix(0,1,(neq+1))
for(i in 1:(neq+1)){
tot11 = 0
tempo1=tempo2=tempo3=0

  for (t in 1:y.power){
    tempo0 <- bjr[t,i]*t*y^(t-1)
    tempo2 <- tempo2+tempo0
  }

if (zy.inter){   for (t in 1:nsoc){
    tempo0 <- hjt[t,i]*Z[,t+1]
    tempo3 <- tempo3+tempo0
  }}
  
if (py.inter){  for (k in 1:(neq+1)){
    tempo0 <- bjk[k,i]*P[,k]
    tempo1 <- tempo1+tempo0
  }}
  
  tot11 <- tot11+tempo1+tempo2+tempo3
  ER[i]=mean(tot11)
}

colnames(ER) <- c(labels.share[1:neq],"others")

## Calculation of standard deviations of income elasticities
## (delta method) if "ERDELTA=TRUE"
if (ERDELTA){
ttt <- colnames(summary(fit3sls)$coefCov)
ER_SE=matrix(0,1,neq+1)
ELASTINCOME_SE=matrix(0,1,neq+1)
for (i in (1:neq)){
tt <- c()
  for (j in (1:y.power)){
 tt <- c(tt,paste("eq",i,"_y",j,sep=""))}


if (zy.inter){
 for (t in (1:nsoc)){
  tt <- c(tt,paste("eq",i,"_yz",t,sep=""))
}}
if (py.inter){
  for (j in 1:neq){
 tt <- c(tt,paste("eq",i,"_ynp",j,sep=""))
}
}

 tnum <- match(tt,ttt)

  DD <- summary(fit3sls)$coefCov[c(tnum),c(tnum)]

 MAT <- c()
 for(r in 1:y.power){
 MAT <- cbind(MAT,r*y^{r-1})}
 if (zy.inter) MAT <- cbind(MAT,Z[,-1])
 if (py.inter) MAT <- cbind(MAT,P[,(1:neq)])

   ER_SE[1,i] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))
ELASTINCOME_SE [1,i] <- ER_SE[1,i]/mean(shares[,i])



}
 ER_SE[1,neq+1] <- -sqrt(sum(ER_SE[1,1:neq]^2))
ELASTINCOME_SE [1,neq+1] <- ER_SE[1,neq+1]/mean(shares[,neq+1])

ER_SE <- as.matrix(ER_SE)
ELASTINCOME_SE <- as.matrix(ELASTINCOME_SE)

colnames(ER_SE) <-colnames(ELASTINCOME_SE) <- c(labels.share[1:neq],"others")
}



}


if (type=="demographics"){
## Calculation of sociodemographic elasticities of budget shares 
## page 13 formula 23 "EASI made EASIER" (Pendakur 2008)
EZ=matrix(0,nsoc,(neq+1))
a <- my.array
for(i in 1:(neq+1)){
tempo4=tempoo=0
  for (t in 1:nsoc){

if (interact){
for (k in 1:(neq+1)){
  if (t %in% interpz){
tempoo <- a[t,k,i]*P[,k]
  tempo4 <- tempo4+tempoo   
}
 }
}   
    tot12 <- gjt[t,i]+hjt[t,i]*y

    tot12 <- tot12+tempo4
    EZ[t,i]=mean(tot12)
  }
}

colnames(EZ) <- c(labels.share[1:neq],"others")
rownames(EZ) <- labels.soc

## Calculation of standard deviations of sociodemographic elasticities 
## (delta method) if "EZDELTA=TRUE"
if (EZDELTA){
ttt <- colnames(summary(fit3sls)$coefCov)
EZ_SE=matrix(0,nsoc,neq+1)
for (i in 1:neq){
  for (j in 1:nsoc){
tt <- c()
    tt <- c(tt,paste("eq",i,"_z",j,sep=""))
if (zy.inter){
  tt <- c(tt,paste("eq",i,"_yz",j,sep=""))
}
if (pz.inter){
      if (j  %in% interpz){
        for (t in 1:neq){
  tt <- c(tt,paste("eq",i,"_np",t,"z",j,sep=""))
}}}


 tnum <- match(tt,ttt)

  DD <- summary(fit3sls)$coefCov[tnum,tnum]

  MAT <- Z[,1]
 
  if (zy.inter)
 MAT <- cbind(MAT,y)
 if (pz.inter){
   if (j  %in% interpz){
 MAT <- cbind(MAT,P[,(1:neq)])
}}
   

   EZ_SE[j,i] <- median(sqrt(diag(as.matrix(MAT)%*%as.matrix(DD)%*%t(as.matrix(MAT)))))

}}
for (j in 1:nsoc) EZ_SE[j,neq+1] <- sqrt(sum(EZ_SE[j,1:neq]^2))

colnames(EZ_SE) <- c(labels.share[1:neq],"others")
rownames(EZ_SE) <- labels.soc


}

}

if (type=="demographics") ER=ELASTINCOME=EP=EPS=EPQ=ELASTPRICE="Not calculated"
if (type=="income") EZ=EP=EPS=EPQ=ELASTPRICE="Not calculated"
if (type=="price") ER=ELASTINCOME=EZ="Not calculated"
if (!EPDELTA)EP_SE=ELASTPRICE_SE="Not calculated"
if (!ERDELTA) ER_SE=ELASTINCOME_SE="Not calculated"
if (!EZDELTA) EZ_SE="Not calculated"


Result <-  list(
#EL1=`*** Semi-elasticities of budget shares in respect to log.prices ***`
EP=EP,
#EL2="*** Standard deviations of Semi-elasticities of budget shares in respect to log.prices ***",
EP_SE=EP_SE,
#EL3="*** Matrix of compensated quantity derivatives with respect to unlogged log.prices ***",
EPS=EPS,
#EL4="*** Compensated (good-specific) expenditures with respect to log.prices ***"
EPQ=EPQ,
#EL5="*** Elasticities of quantities in respect to log.prices ***",
ELASTPRICE=ELASTPRICE,
#EL6="*** Elasticities of quantities in respect to income ***",
ELASTINCOME=ELASTINCOME,
#EL7="*** Standard deviations of Elasticities of quantities in respect to log.prices ***",
ELASTPRICE_SE=ELASTPRICE_SE,
#EL8="*** Standard deviations of Elasticities of quantities in respect to income ***",
ELASTINCOME_SE=ELASTINCOME_SE,
#EL9="*** Semi-elasticities of budget shares in respect to real expenditures ***",
ER=ER,
#EL10="*** Standard deviations of Semi-elasticities of budget shares in respect to real expenditures ***",
ER_SE=ER_SE,
#EL11="*** Semi-elasticities of budget shares in respect to demographics ***",
EZ=EZ,
#EL12="*** Standard deviations of Semi-elasticities of budget shares in respect to demographics ***",
EZ_SE=EZ_SE
)

Result

}
