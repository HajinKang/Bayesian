# Explaining the Gibbs Sampler
# by GEORGE CASELLA and EDWARD I. GEORGE
library(LearnBayes)


# ### 1, 3

library(VGAM)
library(ggplot2)
set.seed(2019121102)
# set.seed(123)
BetaBinomialDraws <- rbetabinom.ab(500, 16, 2, 4, .dontuse.prob = NULL)
ggplot(as.data.frame(BetaBinomialDraws), aes(x=BetaBinomialDraws))+
  geom_histogram()

m <- 500
n <- 16
a <- 2 
b <- 4 
k <- 10
BetaBinomialDrawsGibbs <- numeric(m)
for(i in 1:m){
  x <- 0 
  y <- rbeta(1, 1, 1)
  for(j in 1:k){
    x <- rbinom(1, n, y)
    y <- rbeta(1, shape1 = x+a,shape2 =  b+n-x)
  }
  BetaBinomialDrawsGibbs[i] <- x
}
ggplot(as.data.frame(BetaBinomialDrawsGibbs), aes(x=BetaBinomialDrawsGibbs))+
  geom_histogram()

install.packages("plotrix")
library(plotrix)
l <- list(BetaBinomialDrawsGibbs, BetaBinomialDraws)
multhist(l, breaks = seq(-0.5,16.5, by=1), ylim=c(0,70))
legend("topright", c("Gibbs", "Directily Draws"), col=c("Black", "Gray"), pch=15)







# ### 2, 4
set.seed(123)

X = rep(0, 500)
Y = rep(0, 500)

k = 15
B = 5
for (i in 1:500) {
  x = rep(1, k)
  y = rep(1, k)
  
  for (j in 2:k) {
    temp_x = B+1
    while(temp_x>B) {
      x[j] = rexp(1,y[j-1])
      temp_x = x[j]
    }
    
    temp_y = B+1
    while(temp_y>B) {
      y[j] = rexp(1,x[j])
      temp_y = y[j]
    }     
  }
  
  X[i] = x[k]
  Y[i] = y[k] 
}
print(max(X))
print(max(Y))

hist(X, breaks=100, freq=F)





library(SMPracticals)

# exp.gibbs(u1 = NULL, u2 = NULL, B, I = 15, S = 500)
set.seed(2019121102)
add.exp.lines <- function( exp.out, i, B=5)
{
  dexp.trunc <- function( u, lambda, B ) 
    dexp(u, rate=lambda)/(1-exp(-lambda*B))
  S <- dim(exp.out)[2]
  I <- dim(exp.out)[3]
  u <- seq(0.0001,B,length=1000)
  fu <- rep(0,1000)
  for (s in 1:S) fu <- fu + dexp.trunc(u,exp.out[3-i,s,I],B)/S
  lines(u,fu,col="Pink")
  invisible()
}
par(mfrow=c(2,1))
B <-5; I <- 15; S <- 500
exp.out <- exp.gibbs(B=B,I=I,S=S)
hist(exp.out[1,,I],prob=TRUE,nclass=50,xlab="u1",ylab="PDF",xlim=c(0,2.5),ylim=c(0,0.2))
add.exp.lines(exp.out,1,B=5)
hist(exp.out[2,,I],prob=TRUE,nclass=50,xlab="u2",ylab="PDF",xlim=c(0,B),ylim=c(0,0.2))
add.exp.lines(exp.out,2,B=5)


hist(exp.out[1,,I],prob=FALSE,nclass=50,xlab="x",ylab="PDF",xlim=c(0,B), main="Histogram of X")
add.exp.lines(exp.out,1)
hist(exp.out[1,,I],prob=TRUE,breaks=50,xlab="x",ylab="PDF",xlim=c(0,B),ylim=c(0,1), main="Histogram of X")
add.exp.lines(exp.out,1)


h <- hist(exp.out[1,,I], breaks = 100, plot=FALSE)
h$counts=h$counts/sum(h$counts)
xseq = seq(0,2.5,0.4)
yseq = seq(0,0.16,0.02)
plot(h)
lines(density(exp.out[1,,I]))


X<-exp.out[1,,I]
my.hist <- hist(X, breaks=50, xaxt='n',yaxt='n', ylab="Prob",xlim = c(0,2.5))
ticks <- seq(par("yaxp")[1], par("yaxp")[2], length.out=par("yaxp")[3]+1)
l <- length(X)
max.prob <- max(my.hist$counts)/l
# tick.labels <- head(pretty(c(0, max.prob)), -1)
tick.labels <- pretty(c(0, max.prob))
ticks <- tick.labels * l
print(tick.labels)
print(sum(my.hist$counts/l))
axis(side=1,at=seq(0,2.5,0.4))
axis(side=2, at=ticks, labels=tick.labels)


xi = seq(0.001, 5, length.out=50)
B = 5
Tmarginal=(1-exp(-B*xi))/xi
TTmarginal=Tmarginal/sum(Tmarginal)
lines(xi,TTmarginal, col="Red", lwd=2, lty=2)


# 
fhat<-numeric(I)
for (i in 1:I){
  fhat[i]<-mean(exp.out[1,,i])
}
order(fhat)
Tfhat<-fhat/sum(fhat)
# 