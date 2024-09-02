knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
```{r code=xfun::read_utf8('../Pascal2021/dr09.out'), echo=TRUE, eval=FALSE}
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
### Special implementations
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
```{r code=xfun::read_utf8('../Pascal2021/dr03.out'), echo=TRUE, eval=FALSE}
$$  R_n x = z_1 $$
and the minimal sum of squares turns out to be the cross-product $z_2' z_2$.
library('nlsr')
citation('nslr')
citation('nlsr')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr09.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr09.R')
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
```{r code=xfun::read_utf8('../Pascal2021/dr03.out'), echo=TRUE, eval=FALSE}
```{r code=xfun::read_utf8('../Pascal2021/dr0506.out'), echo=TRUE, eval=FALSE}
getwd()
cubfn <- function(x) { x*(x*x-1)-5}
curve(cubfn)
res <- optim(par=0.0, fn=cubfn, method="Brent", lower=c(0), upper=c(1))
res
?curve
curve(cubfn, from=-10, to=10)
curve(cubfn, from=-5, to=5)
curve(cubfn, from=-3, to=3)
curve(cubfn, from=-2, to=2)
curve(cubfn, from=-1.5, to=1.5)
cubfn(1)
cubfn(2.5)
cubfn(2.5)
cubfn <- function(x) { x*(x*x-2)-5}
curve(cubfn, from=-1.5, to=1.5)
res <- optim(par=0.0, fn=cubfn, method="Brent", lower=c(0), upper=c(1))
res
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
fminfn <-function(x){
val<-(x[2]-x[1]^2)^2)*100.0+(1.0-x[1])^2
fminfn <-function(x){
val<-((x[2]-x[1]^2)^2)*100.0+(1.0-x[1])^2
val
}
x0<-c(-1.2,1)
cat("Check fn at c(-1.2,1)=",fminfn(x0),"\n")
rslt <- hjn(par=x0, fn=fminfn, control=list(trace=1))
require("../R/hjn.R")
fminfn <-function(x){
val<-((x[2]-x[1]^2)^2)*100.0+(1.0-x[1])^2
val
}
x0<-c(-1.2,1)
cat("Check fn at c(-1.2,1)=",fminfn(x0),"\n")
rslt <- hjn(par=x0, fn=fminfn, control=list(trace=1))
library("../R/hjn.R")
source("../R/hjn.R")
fminfn <-function(x){
val<-((x[2]-x[1]^2)^2)*100.0+(1.0-x[1])^2
val
}
x0<-c(-1.2,1)
cat("Check fn at c(-1.2,1)=",fminfn(x0),"\n")
rslt <- hjn(par=x0, fn=fminfn, control=list(trace=1))
print(rslt)
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr03.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr03.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr03.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr03ca.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/dr03ca.R')
d <- c(5.0489173, 0.6431041, 0.3079785)
V <- matrix(c( 0.3279853, -0.7369762, 0.591009, 0.591009, -0.3279853, -0.7369762,
0.7369762, 0.591009, 0.3279853), nrow=3, ncol=3)
V
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/zz-unversioned/testbasa14.R')
A
A %*% V
Dm <- diag(d)
Dm
AV <- A %*% V
Dm <- diag(d)
Dmv <- Dm %*% V
AV - Dmv
d <- c(5.0489173, 0.6431041, 0.3079785)
V <- matrix(c( 0.3279853, -0.7369762, 0.591009, 0.591009, -0.3279853, -0.7369762,
0.7369762, 0.591009, 0.3279853), nrow=3, ncol=3)
A <- matrix(0.0, nrow=3, ncol=3)
for (i in 1:3){
for (j in 1:3){
A[i,j] <- min(i,j)-2.0
}
if (i <= j) {A[i,i] <- i}
}
Vsave <- V
V <- t(V)
AV <- A %*% V
Dm <- diag(d)
Dmv <- Dm %*% V
AV - Dmv
eig(A)
eigen(A)
V
reig <- eigen(A)
Vr <- reig$vectors
dr <- reig$values
AVr <- A %*% Vr
Dmr <- diag(dr)
Dmrv <- Dmr %*% Vr
AVr - Dmrv
reig <- eigen(A)
Vr <- reig$vectors
dr <- reig$values
AVr <- A %*% Vr
Dmr <- diag(dr)
Dmrv <- Vr %*% Dmr
AVr - Dmrv
AV <- A %*% V
Dm <- diag(d)
Dmv <- V %*% Dm
AV - Dmv
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
?system.time
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/Rtimetest.R')
source('/j19z/j19store/versioned/Nash-Compact-Numerical-Methods/R/Rtimetest.R')
knitr::opts_chunk$set(echo = TRUE)
## ??require(bookdown) # language engine to display text??
install.packages("rticles")
library(rticles)
install.packages("distill")
library(distill)
library(pracma)
y <- c(20.91676, 20.65219, 20.39272, 20.58692, 21.64712, 23.30965, 23.35657,
24.22724, 24.83439, 24.34865, 23.13173, 21.96117)
t <- c(1, 2, 3, 4 , 5 , 6, 7, 8, 9, 10, 11, 12)
# Fitting function
fit <- function(x, a, b, c) {a+b*sin(2*pi*x)+c*cos(2*pi*x)}
res <- nls(y ~ fit(t, a, b, c), data=data.frame(t,y), start = list(a=1,b=0,
c=1))
f1 <- y ~  a+b*sin(2*pi*x)+c*cos(2*pi*x)
res1 <- nls(f1, data=data.frame(t,y), start=list(a=1,b=0, c=1))
f1 <- y ~  a+b*sin(2*pi*t)+c*cos(2*pi*t)
res1 <- nls(f1, data=data.frame(t,y), start=list(a=1,b=0, c=1))
library(nlsr)
res1n <- nlxb(f1, data=data.frame(t,y), start=list(a=1,b=0, c=1))
res1n
rm(list=ls())
ls()
?pracma
?wilkindon
?wilkinson
wilkinson(7)
wilkinson
n<-8
m<-(n-1)/2
m
Diag(abs(-m:m))
Diag(3,1)
wilkinson(8)
source('~/temp/wminus.R')
wminus(7)
A<-wminus(100)
A<-A( ,1:50)
A<-A[ ,1:50]
?QR
??QR
myqr<-qr(A)
myqr
str(myqr)
diag(qr$qr)
diag(myqr$qr)
dd<-diag(myqr$qr)
?product
??product
?cumprod(dd)
cumprod(dd)
source('~/temp/Givenswork2103/Give210328.R', echo=TRUE)
mdata<-c(1.00,   3.00,  5.00,
1.10,   2.90,  10.00,
1.20,   2.80,   8.0,
1.30,   2.70,  12.0,
1.40,   2.60,   1.0)
mdata
A <-matrix(A, ncol=3)
A <-matrix(mdata, ncol=3)
svd(A)
tdat<-c(2.7018512172212588, 6.1809473051500028, 15.766967377208989,
-0.0000000000000000, 1.0468478451804264, 3.7686522426495328,
0.0000000000000000, 9.4001039556654356E-018,   8.4380092438915941 )
B<-matrix(tdat,ncol=3)
svd(B)
savehistory("210331A.txt")
aqr<-qr(A)
aqr
?qr
?qr.qr
??qr
?qr.R
qr.R(A)
qr.R(qra)
ls()
qr.R(aqr)
raqr<-qr.R(aqr)
svd(raqr)
savehistory("210331A.txt")
ls()
svd(A)$d
svd(qr.R(qr(A)))$d
A
getwd()
rm(list=ls())
mdata<-c(1.00,   3.00,  5.00,
1.10,   2.90,  10.00,
1.20,   2.80,   8.0,
1.30,   2.70,  12.0,
1.40,   2.60,   1.0)
mdata
A <-matrix(A, ncol=3)
A <-matrix(mdata, ncol=3, byrow=TRUE)
svd(A)
A <-matrix(mdata, ncol=3, byrow=TRUE)
A
svd(A)
tdat<-c(2.7018512172212588, 6.1809473051500028, 15.766967377208989,
-0.0000000000000000, 1.0468478451804264, 3.7686522426495328,
0.0000000000000000, 9.4001039556654356E-018,   8.4380092438915941 )
B<-matrix(tdat,ncol=3, byrow=TRUE)
B
svd(B)
savehistory("210331A.txt")
nsvdres<-c(19.267636745626650, 3.0317062066817919, 0.40857277802700392)
tsvdres<-svd(B)
tsvdres-nsvdres
str(tsvdres)
tsvdres<-svd(B)$d
tsvdres-nsvdres
ls()
asvdres<-svd(A)$d
asvdres-nsvdres
dsvdres<-c(19.267636745626646, 3.0317062066817910,  0.40857277802700365)
dsvdres-asvdres
dsvdres-tstvdres
dsvdres-tsvdres
dsvdres-nsvdres
savehistory("210331B.R")
