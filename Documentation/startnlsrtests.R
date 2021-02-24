# set parameters
a <- 1
b <- 2
c <- 3
np <- 15
tt <-1:np
yy <- 100*a/(1+10*b*exp(-0.1*c*tt))
plot.new()
plot(tt, yy, type='l')
set.seed(123456)
ev <- runif(np)
ev <- ev - mean(ev)
y1 <- yy + ev
points(tt,y1,type='p', col="green")
y2 <- yy + 5*ev
points(tt,y2,type='p', col="blue")
y3 <- yy + 10*ev
lg3d15 <- data.frame(tt, yy, y1, y2, y3)
points(tt,y3,type='p', col="red")
library(nlsr)
sol0 <- nlxb(yy ~ a0/(1+b0*exp(-c0*tt)), data=lg3d15, start=list(a0=1, b0=1, c0=1))
print(sol0)
sol1 <- nlxb(y1 ~ a1/(1+b1*exp(-c1*tt)), data=lg3d15, start=list(a1=1, b1=1, c1=1))
print(sol1)
sol2 <- nlxb(y2 ~ a2/(2+b2*exp(-c2*tt)), data=lg3d15, start=list(a2=1, b2=1, c2=1))
print(sol2)
sol3 <- nlxb(y3 ~ a3/(3+b3*exp(-c3*tt)), data=lg3d15, start=list(a3=1, b3=1, c3=1))
print(sol3)
# yy <- 100*a/(1+10*b*exp(-0.1*c*tt))
# plot(tt, yy, type='l')
# set.seed(123456)
# ev <- runif(np)
# ev <- ev - mean(ev)
# y1 <- yy + ev
# points(tt,y1,type='p', col="green")
# y2 <- yy + 5*ev
# points(tt,y2,type='p', col="blue")
# y3 <- yy + 10*ev
# lg3d15 <- data.frame(tt, yy, y1, y2, y3)
# points(tt,y3,type='p', col="red")
# library(nlsr)
# sol0 <- nlxb(yy ~ a0/(1+b0*exp(-c0*tt)), data=lg3d15)
# print(sol0)
# sol1 <- nlxb(y1 ~ a1/(1+b1*exp(-c1*tt)), data=lg3d15)
# print(sol1)
# sol2 <- nlxb(y2 ~ a2/(2+b2*exp(-c2*tt)), data=lg3d15)
# print(sol2)
# sol3 <- nlxb(y3 ~ a3/(3+b3*exp(-c3*tt)), data=lg3d15)
# print(sol3)
## want to change number of points
# np <- 150
# tt <- (1:np)/10
# yy <- 100*a/(1+10*b*exp(-0.1*c*tt))
# set.seed(123456)
# ev <- runif(np)
# ev <- ev - mean(ev)
# y1 <- yy + ev
# y2 <- yy + 5*ev
# y3 <- yy + 10*ev
# lg3d150 <- data.frame(tt, yy, y1, y2, y3)
# np <- 1500
# tt <- (1:np)/100
# yy <- 100*a/(1+10*b*exp(-0.1*c*tt))
# set.seed(123456)
# ev <- runif(np)
# ev <- ev - mean(ev)
# y1 <- yy + ev
# y2 <- yy + 5*ev
# y3 <- yy + 10*ev
# lg3d1500 <- data.frame(tt, yy, y1, y2, y3)
# f0 <- yy ~ a0/(1+b0*exp(-c0*tt))
# f1 <- y1 ~ a1/(1+b1*exp(-c1*tt))
# f2 <- y2 ~ a2/(2+b2*exp(-c2*tt))
# f3 <- y3 ~ a3/(3+b3*exp(-c3*tt))
