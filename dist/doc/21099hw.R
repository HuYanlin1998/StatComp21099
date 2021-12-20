## -----------------------------------------------------------------------------
library(ggplot2)
iris.summary <- aggregate(iris[1:4], list(iris$Species), mean)
names(iris.summary)[1] <- 'Species'
ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, col = Species)) + geom_point() + geom_point(data = iris.summary, shape = 15, size = 5)

## -----------------------------------------------------------------------------
head(iris)
library(kableExtra)
data<-head(iris)
data %>%
  kbl(caption = "table of iris") %>%
  kable_classic(full_width = F, html_font = "Cambria")

## -----------------------------------------------------------------------------
n <- 2000
u <- runif(n)
sigma=1   # sigma=1
x=sigma*sqrt(-2*log(1-u))
hist(x, prob = TRUE, breaks=c(0:50)/10, main = expression(f(x)==x*e^{-x^2/2}))
y <- seq(0, 5, .01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
n <- 2000
u <- runif(n)
sigma=2   # sigma=2
x=sigma*sqrt(-2*log(1-u))
ran<-floor(range(x)[2])+1
hist(x, prob = TRUE, breaks=c(0:(ran*10))/10, main = expression(f(x)==x/4*e^{-x^2/8}))
y <- seq(0, ran, .01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))

## -----------------------------------------------------------------------------
n <- 2000
u <- runif(n)
sigma=0.5   # sigma=0.5
x=sigma*sqrt(-2*log(1-u))
ran<-floor(range(x)[2])+1
hist(x, prob = TRUE, breaks=c(0:(ran*10))/10, main = expression(f(x)==4*x*e^{-2*x^2}))
y <- seq(0, ran, .01)
lines(y,(y/sigma^2)*exp(-y^2/(2*sigma^2)))


## -----------------------------------------------------------------------------
p1<-0.75
p2<-1-p1
n<-1000
u<-runif(n)
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
X<-rep(0,n)
for(i in 1:n){
  if(u[i]<=p1){
    X[i]<-X1[i]
  }
  else{X[i]<-X2[i]}
}
ran1<-floor(range(X)[1])
ran2<-floor(range(X)[2])+1
hist(X,prob = TRUE, breaks=c((ran1*10):(ran2*10))/10,main = expression(p1==0.75))
y <- seq(ran1, ran2, .01)
lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1))


## -----------------------------------------------------------------------------
p1<-0.5
p2<-1-p1
n<-1000
u<-runif(n)
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
X<-rep(0,n)
for(i in 1:n){
  if(u[i]<=p1){
    X[i]<-X1[i]
  }
  else{X[i]<-X2[i]}
}
ran1<-floor(range(X)[1])
ran2<-floor(range(X)[2])+1
hist(X,prob = TRUE, breaks=c((ran1*10):(ran2*10))/10,main = expression(p1==0.5))
y <- seq(ran1, ran2, .01)
lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1))


## -----------------------------------------------------------------------------
p1<-0.25
p2<-1-p1
n<-1000
u<-runif(n)
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
X<-rep(0,n)
for(i in 1:n){
  if(u[i]<=p1){
    X[i]<-X1[i]
  }
  else{X[i]<-X2[i]}
}
ran1<-floor(range(X)[1])
ran2<-floor(range(X)[2])+1
hist(X,prob = TRUE, breaks=c((ran1*10):(ran2*10))/10,main = expression(p1==0.25))
y <- seq(ran1, ran2, .01)
lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1))


## -----------------------------------------------------------------------------
p1<-0.1
p2<-1-p1
n<-1000
u<-runif(n)
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
X<-rep(0,n)
for(i in 1:n){
  if(u[i]<=p1){
    X[i]<-X1[i]
  }
  else{X[i]<-X2[i]}
}
ran1<-floor(range(X)[1])
ran2<-floor(range(X)[2])+1
hist(X,prob = TRUE, breaks=c((ran1*10):(ran2*10))/10,main = expression(p1==0.1))
y <- seq(ran1, ran2, .01)
lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1))


## -----------------------------------------------------------------------------
p1<-0.9
p2<-1-p1
n<-1000
u<-runif(n)
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
X<-rep(0,n)
for(i in 1:n){
  if(u[i]<=p1){
    X[i]<-X1[i]
  }
  else{X[i]<-X2[i]}
}
ran1<-floor(range(X)[1])
ran2<-floor(range(X)[2])+1
hist(X,prob = TRUE, breaks=c((ran1*10):(ran2*10))/10,main = expression(p1==0.9))
y <- seq(ran1, ran2, .01)
lines(y,p1*dnorm(y,0,1)+p2*dnorm(y,3,1))


## -----------------------------------------------------------------------------
gam<-function(lam=1,alfa=1,gamm=1,t=10,n=1000){
  #lam is poisson(lam) parameter,alfa is shape parameter,gamm is scale parameter,n is sample size,
  #t is time
  Xt<-rep(0,n)
  Nt<-rpois(n,lam*t)
  for(i in 1:n){
    if(Nt[i]==0){
      Xt[i]<-0
    }
    else{Xt[i]<-sum(rgamma(Nt[i],shape=alfa,scale=gamm))}
  }
  bar<-mean(Xt)
  fangcha<-var(Xt)
  bar_the<-lam*t*alfa*gamm
  fangcha_the<-lam*t*(alfa^2+alfa)*gamm^2
  cat("sample mean at time t=",t,"is",bar,"\n")
  cat("sample variance at time t=",t,"is",fangcha,"\n")
  cat("real mean at time t=",t,"is",bar_the,"\n")
  cat("real variance at time t=",t,"is",fangcha_the,"\n")
}  
gam(1,1,1,10,1000)
gam(2,1,1,10,1000)
gam(1,2,3,10,1000)

## -----------------------------------------------------------------------------
pbeta1<-function(x=0.5,a=3,b=3,n=1000){
  #a,b,are two shape parameters,n is sample size, x is in (0,1),this function return the value of cdf of beta(a,b) at x
  beta<-beta(a,b)
  u<-runif(n,0,x)
  Fx<-(x/beta)*mean(u^(a-1)*(1-u)^(b-1))
  return(Fx)
}
# next we Compare the estimates with the values returned by the pbeta function in R
abs(pbeta1(0.1,3,3,3000)-pbeta(0.1,3,3))
abs(pbeta1(0.2,3,3,3000)-pbeta(0.2,3,3))
abs(pbeta1(0.3,3,3,3000)-pbeta(0.3,3,3))
abs(pbeta1(0.4,3,3,3000)-pbeta(0.4,3,3))
abs(pbeta1(0.5,3,3,3000)-pbeta(0.5,3,3))
abs(pbeta1(0.6,3,3,3000)-pbeta(0.6,3,3))
abs(pbeta1(0.7,3,3,3000)-pbeta(0.7,3,3))
abs(pbeta1(0.8,3,3,3000)-pbeta(0.8,3,3))
abs(pbeta1(0.9,3,3,3000)-pbeta(0.9,3,3))

## -----------------------------------------------------------------------------
rRayleigh<-function(sigma=1,antithetic=TRUE,m=3000){
  u<-runif(m)
  v<-runif(2*m)
  if(antithetic==TRUE){
    x<-(sigma*sqrt(-2*log(1-u))+sigma*sqrt(-2*log(u)))/2
  } else{
    x<-(sigma*sqrt(-2*log(1-v[1:m]))+sigma*sqrt(-2*log(1-v[(m+1):(2*m)])))/2
  }
  return(x)
}


var(rRayleigh())
var(rRayleigh(antithetic=FALSE))

abs(var(rRayleigh(antithetic=FALSE))-var(rRayleigh()))/var(rRayleigh(antithetic=FALSE))

## -----------------------------------------------------------------------------
Important1<-function(n=3000){
  u<-runif(n)
  x<-sqrt(1-2*log(u))
  theta<-mean(x)/sqrt(2*exp(1)*pi)
  return(theta)
}
Important1()


## -----------------------------------------------------------------------------

n<-20
quant<-qt(0.975,n-1)
num<-0
for(i in 1:1000){
  x<-rchisq(n,2)
  m<-mean(x)
  se<-sqrt(var(x))
  if((m-quant*se/sqrt(n)<2)&(m+quant*se/sqrt(n)>2)){
    num<-num+1
  }
}
mean_rate_chisq<-num/1000
mean_rate_chisq

num<-0
for(i in 1:1000){
  x<-rnorm(n,0,2)
  m<-mean(x)
  se<-sqrt(var(x))
  if((m-quant*se/sqrt(n)<0)&(m+quant*se/sqrt(n)>0)){
    num<-num+1
  }
}
mean_rate_norm<-num/1000
mean_rate_norm

num<-0
alpha <- .05
for(i in 1:1000){
  x <- rnorm(n, mean=0, sd=2)
  UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
  if(UCL>4){
    num<-num+1
  }
}
fangcha_rate_norm<-num/1000
fangcha_rate_norm

num<-0
alpha <- .05
for(i in 1:1000){
  x <-rchisq(n,2)
  UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1)
  if(UCL>4){
    num<-num+1
  }
}
fangcha_rate_chisq<-num/1000
fangcha_rate_chisq

## -----------------------------------------------------------------------------
m <- 1e5
n <- 20
alpha<-0.05
u0=1
p<-rep(0,m)
for(i in 1:m){
  x<-rchisq(n,1)
  m<-mean(x)
  se<-sqrt(var(x))
  t0<-sqrt(n)*(m-u0)/se
  p[i]<-2*(1-pt(abs(t0),n-1))
}

mean(p<=0.05)

## -----------------------------------------------------------------------------
m <- 1e5
n <- 20
alpha<-0.05
u0=1
p<-rep(0,m)
for(i in 1:m){
  x<-runif(n,0,2)
  me<-mean(x)
  se<-sqrt(var(x))
  t0<-sqrt(n)*(me-u0)/se
  p[i]<-2*(1-pt(abs(t0),n-1))
}
mean(p<=0.05)

## -----------------------------------------------------------------------------
m <- 1e5
n <- 20
alpha<-0.05
u0=1
p<-rep(0,m)
for(i in 1:m){
  x<-rexp(n,1)
  me<-mean(x)
  se<-sqrt(var(x))
  t0<-sqrt(n)*(me-u0)/se
  p[i]<-2*(1-pt(abs(t0),n-1))
}
mean(p<=0.05)


## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)
#  n <- c(10, 20, 30, 50, 100, 500)
#  d=2                                #degrees of freedom
#  cv2 <- qchisq(.975,d*(d+1)*(d+2)/6) #crit. values
#  cv1 <- qchisq(.025,d*(d+1)*(d+2)/6)
#  
#  sk <- function(x) {
#    #computes the sample skewness coeff.
#    x<-scale(x,center = T,scale=F)
#    sigmahat<-var(x)
#    nn<-nrow(x)
#    b1d<-sum(colSums((x%*%solve(sigmahat)%*%t(x))^3))
#    b1d<-b1d/6/nn
#    return( b1d)
#  }
#  
#  
#  p.reject <- numeric(length(n)) #to store sim. results
#  m <- 500  #num. repl. each sim.
#  for (i in 1:length(n)) {
#    sktests <- numeric(m) #test decisions
#    for (j in 1:m) {
#      x <- mvrnorm(n[i], rep(0, d),diag(d))
#      #test decision is 1 (reject) or 0
#      sktests[j] <- as.integer(sk(x)>= cv2 |sk(x)<= cv1 )
#    }
#    p.reject[i] <- mean(sktests) #proportion rejected
#  }
#  p.reject

## ----eval=FALSE---------------------------------------------------------------
#  library(MASS)
#  alpha <- .1
#  n <- 30
#  m <- 500
#  epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
#  N <- length(epsilon)
#  pwr <- numeric(N)
#  sk <- function(x) {
#    #computes the sample skewness coeff.
#    x<-scale(x,center = T,scale=F)
#    sigmahat<-var(x)
#    nn<-nrow(x)
#    b1d<-sum(colSums((x%*%solve(sigmahat)%*%t(x))^3))
#    b1d<-b1d/6/nn
#    return( b1d)
#  }
#  #critical value for the skewness test
#  d=2                                #degrees of freedom
#  cv2 <- qchisq(1-alpha/2,d*(d+1)*(d+2)/6) #crit. values
#  cv1 <- qchisq(alpha/2,d*(d+1)*(d+2)/6)
#  #generating mixed multivariate normal
#  hunhe<-function(size=30,dimension=2,eps=0.1){
#    xx<-matrix(rep(0,size*dimension),size,dimension)
#    comp <- sample(c(0,1), replace = TRUE, size = n, prob = c(1-eps, eps))
#    for(k in 1:size){
#      if(comp[k]==0){
#        xx[k,]<-mvrnorm(1, rep(0,d), diag(d))
#      }
#      else {xx[k,]<-mvrnorm(1, rep(0,d), 10*diag(d))}
#    }
#    return(xx)
#  }
#  
#  for (j in 1:N) { #for each epsilon
#    e <- epsilon[j]
#    sktests <- numeric(m)
#    for (i in 1:m) { #for each replicate
#      x<-hunhe(n,d,e)
#      sktests[i] <- as.integer(sk(x) >= cv2|sk(x) <= cv1)
#    }
#    pwr[j] <- mean(sktests)
#  }
#  #plot power vs epsilon
#  
#  plot(epsilon, pwr, type = "b",
#       xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = .1, lty = 3)
#  se <- sqrt(pwr * (1-pwr) / m) #add standard errors
#  lines(epsilon, pwr+se, lty = 3)
#  lines(epsilon, pwr-se, lty = 3)
#  

## -----------------------------------------------------------------------------

library(bootstrap)
R<-1000
n<-nrow(scor)
theta_hat_b<-rep(0,R)
theta_hat<-(eigen((n-1)/n*var(scor))$values[1])/sum(eigen((n-1)/n*var(scor))$values)
for(i in 1:R){
  ind<-sample(c(1:n),replace = TRUE)
  theta<-eigen((n-1)/n*var(scor[ind,]))
  theta_hat_b[i]<-(theta$values[1])/(sum(theta$values))
}
theta_hat_bias<-mean(theta_hat_b)-theta_hat
theta_hat_bias
theta_hat_se<-sd(theta_hat_b)
theta_hat_se

## -----------------------------------------------------------------------------
theta_hat_jack<-rep(0,n)
for(i in 1:n){
  theta_hat_jack[i]<-(eigen((n-1)/n*var(scor[-i,]))$values[1])/sum(eigen((n-1)/n*var(scor[-i,]))$values)
}
theta_hat_bias_jack<-(n-1)*(mean(theta_hat_jack)-theta_hat)
theta_hat_bias_jack
theta_hat_se_jack<-sqrt((n-1)^2/n*var(theta_hat_jack))
theta_hat_se_jack

## -----------------------------------------------------------------------------
alpha<-0.05
low1<-sort(theta_hat_b)[R*alpha/2]
up1<-sort(theta_hat_b)[R-R*alpha/2]
z0<-qnorm(sum(theta_hat_b<theta_hat)/R)
ahat<-sum((theta_hat_b-theta_hat)^3)/6/(sum((theta_hat_b-theta_hat)^2))^(3/2)
alpha1<-pnorm(z0+(z0+qnorm(alpha/2))/(1-ahat*(qnorm(alpha/2))))
alpha2<-pnorm(z0+(z0+qnorm(1-alpha/2))/(1-ahat*(qnorm(1-alpha/2))))
low2<-sort(theta_hat_b)[R*alpha1]
up2<-sort(theta_hat_b)[R*alpha2]
cat("percentile=[",low1,",",up1,"]","\t",
         "BCa=[",low2,",",up2,"]","\t")



## -----------------------------------------------------------------------------
# normal
alpha=0.05
R<-100 # number of resample
N<-100 # num of MC
n<-5  # num of sample
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}
low_norm<-up_norm<-rep(0,N)
low_basic<-up_basic<-rep(0,N)
low_percent<-up_percent<-rep(0,N)


for(i in 1:N){
  x<-rnorm(n)
  sk_b<-rep(0,R)
  for(j in 1:R){
    ind<-sample(c(1:n),replace = TRUE)
    sk_b[j]<-sk(x[ind])
  }
  sk_hat<-sk(x)
  sk_hat_se<-sd(sk_b)
  low_norm[i]<-sk_hat-qnorm(1-alpha/2)*sk_hat_se
  up_norm[i]<-sk_hat-qnorm(alpha/2)*sk_hat_se
  low_basic[i]<-2*sk_hat-sort(sk_b)[R-R*alpha/2]
  up_basic[i]<-2*sk_hat-sort(sk_b)[R*alpha/2]
  low_percent[i]<-sort(sk_b)[R*alpha/2]
  up_percent[i]<-sort(sk_b)[R-R*alpha/2]
}
por_left_norm<-sum(low_norm>0)/R
por_right_norm<-sum(up_norm<0)/R
por_norm<-1-por_left_norm-por_right_norm
por_left_basic<-sum(low_basic>0)/R
por_right_basic<-sum(up_basic<0)/R
por_basic<-1-por_left_basic-por_right_basic
por_left_percent<-sum(low_percent>0)/R
por_right_percent<-sum(up_percent<0)/R
por_percent<-1-por_left_percent-por_right_percent
cat("coverage rate of normal BCI=",por_norm,",rate of missing on left=",
    por_left_norm,",rate of missing on right=",por_right_norm,"\n",
    "coverage rate of basic BCI=",por_basic,",rate of missing on left=",
    por_left_basic,",rate of missing on right=",por_right_basic,"\n",
    "coverage rate of percentile BCI=",por_percent,",rate of missing on left=",
    por_left_percent,",rate of missing on the right=",por_right_percent,"\n")


# chisq(5)

for(i in 1:N){
  x<-rchisq(n,5)
  sk_b<-rep(0,R)
  for(j in 1:R){
    ind<-sample(c(1:n),replace = TRUE)
    sk_b[j]<-sk(x[ind])
  }
  sk_hat<-sk(x)
  sk_hat_se<-sd(sk_b)
  low_norm[i]<-sk_hat-qnorm(1-alpha/2)*sk_hat_se
  up_norm[i]<-sk_hat-qnorm(alpha/2)*sk_hat_se
  low_basic[i]<-2*sk_hat-sort(sk_b)[R-R*alpha/2]
  up_basic[i]<-2*sk_hat-sort(sk_b)[R*alpha/2]
  low_percent[i]<-sort(sk_b)[R*alpha/2]
  up_percent[i]<-sort(sk_b)[R-R*alpha/2]
}
skreal<-sqrt(8/5)
por_left_norm<-sum(low_norm>skreal)/R
por_right_norm<-sum(up_norm<skreal)/R
por_norm<-1-por_left_norm-por_right_norm
por_left_basic<-sum(low_basic>skreal)/R
por_right_basic<-sum(up_basic<skreal)/R
por_basic<-1-por_left_basic-por_right_basic
por_left_percent<-sum(low_percent>skreal)/R
por_right_percent<-sum(up_percent<skreal)/R
por_percent<-1-por_left_percent-por_right_percent
cat("coverage rate of normal BCI=",por_norm,",rate of missing on left=",
    por_left_norm,",rate of missing on right=",por_right_norm,"\n",
    "coverage rate of basic BCI=",por_basic,",rate of missing on left=",
    por_left_basic,",rate of missing on right=",por_right_basic,"\n",
    "coverage rate of percentile BCI=",por_percent,",rate of missing on left=",
    por_left_percent,",rate of missing on the right=",por_right_percent,"\n")



## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1)
#  x<-rnorm(10,sd=1)
#  y<-rnorm(10,sd=2)
#  R <- 999    #number of replicates
#  z <- c(x, y)  #pooled sample
#  K <- 1:20
#  reps <- numeric(R) #storage for replicates
#  t0 <- cor(x, y, method = "spearman")
#  for (i in 1:R) {
#    #generate indices k for the first sample
#    k <- sample(K, size =10, replace = FALSE)
#    x1 <- z[k]
#    y1 <- z[-k] #complement of x1
#    reps[i] <- cor(x1, y1,method = "spearman")
#  }
#  
#  phat <- mean(c(t0, reps) >= t0)
#  cat("achieved significance level of the permutation test is " ,phat,"\n")
#  test<-cor.test(x,y,method = "spearman")
#  p<-test$p.value
#  cat("p-value reported by cor.test is ",p,"\n")
#  

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  library(RANN)
#  library(MASS)
#  library(boot)
#  library(energy)
#  library(Ball)
#  
#  # (1) normal with mu1=mu2=0,sd1=1.5,sd2=1
#  
#  
#  Tn <- function(z, ix, sizes,k) {
#    n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#    if(is.vector(z)) z <- data.frame(z);
#    z <- z[ix, ];
#    NN <- nn2(data=z, k=k+1) # What is the first column?
#    block1 <- NN$nn.idx[1:n1,-1]
#    block2 <- NN$nn.idx[(n1+1):n,-1]
#    i1 <- sum(block1 <= n1); i2 <- sum(block2 > n1)
#    (i1 + i2) / (k * n)
#  }
#  
#  
#  
#  m <- 1e3; k<-3; p<-2; mu <- 0;
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)}
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
#    y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,
#                             seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  # (2)normal with mu1=0,mu2=0.3,sd1=1.5,sd2=1
#  
#  m <- 1e3; k<-3; p<-2; mu <- 0.3;
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)}
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
#    y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,
#                             seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  # (3) t distribution with df=1 vs 0.4N(0,1)+0.6N(1,1)
#  
#  m <- 1e3; k<-3; p<-1;
#  n1 <- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)
#  hunhe<-function(samsize,mu1=0,mu2=1,sd1=1,sd2=1,p=0.4){
#    x<-rep(0,samsize)
#    u<-runif(samsize)
#    for(i in 1:samsize){
#      if(u[i]<=0.4){
#        x[i]<-rnorm(1,mu1,sd1)
#      }
#      if(u[i]>0.4){
#        x[i]<-rnorm(1,mu2,sd2)
#      }
#    }
#    return(x)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rt(n1,df=n1),ncol=1);
#    y <- matrix(hunhe(n2),ncol=1);
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,
#                             seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  # (4) Unbalanced samples
#  
#  m <- 1e3; k<-3; p<-2; mu <- 0.3;
#  n1 <-10; n2 <- 100; R<-999; n <- n1+n2; N = c(n1,n2)
#  eqdist.nn <- function(z,sizes,k){
#    boot.obj <- boot(data=z,statistic=Tn,R=R,
#                     sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#    p.value <- mean(ts>=ts[1])
#    list(statistic=ts[1],p.value=p.value)}
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#    x <- matrix(rnorm(n1*p,0,1.5),ncol=p);
#    y <- cbind(rnorm(n2),rnorm(n2,mean=mu));
#    z <- rbind(x,y)
#    p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#    p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#    p.values[i,3] <- bd.test(x=x,y=y,num.permutations=R,
#                             seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow <- colMeans(p.values<alpha)
#  pow
#  
#  

## -----------------------------------------------------------------------------
set.seed(0)
m <- 10000
x <- numeric(m)
#提议分布用正态分布
x[1] <- rnorm(1)
k <- 0
u <- runif(m)
for (i in 2:m) {
  xt <- x[i-1]
  y <- rnorm(1, mean=xt)
  num <- dcauchy(y) 
  den <- dcauchy(xt) 
  if (u[i] <= num/den) x[i] <- y else {
    x[i] <- xt
    k <- k+1 #y is rejected
  } 
}

x<-x[-(1:1000)]
a <- ppoints(100)
QR <- qcauchy(a) #quantiles of Cauchy
Q <- quantile(x, a)
qqplot(QR, Q, main="",
       xlab="Cauchy Quantiles", ylab="Sample Quantiles")
hist(x, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR,dcauchy(QR))

a <- c(.05, seq(.1, .9, .1), .95)
Q <- qcauchy(a)
Qx<-quantile(x,a)
print(round(cbind(Q,Qx),3))

## -----------------------------------------------------------------------------
set.seed(1)
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample

n<-15
a<-5
b<-5
X[1,]<-c(1,0.3)
for (i in 2:N) {
  x2 <- X[i-1, 2]
  X[i, 1] <- rbinom(1, n, x2)
  x1 <- X[i, 1]
  X[i, 2] <- rbeta(1, x1+a, n-x1+b)
}
x <- X[(burn+1):N, ]

plot(x, main="", cex=.5, xlab=bquote(X[1]),
     ylab=bquote(X[2]), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
set.seed(1)
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #row means
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #upper variance est.
  r.hat <- v.hat / W #G-R statistic
  return(r.hat)
}
cauchy.chain <- function(N, X1) {
  #generates a Metropolis chain for Cauchy(0,1)
  #with N(X[t], 1) proposal distribution
  #and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, mean=xt,sd=1)
    num <- dcauchy(y) 
    den <- dcauchy(xt) 
    if (u[i] <= num/den) x[i] <- y else {
      x[i] <- xt
      k <- k+1 #y is rejected
    } 
  }
  return(x)
}
k <- 4 #number of chains to generate
n <- 20000 #length of chains
b <- 1000 #burn-in length
#choose overdispersed initial values
x0 <- c(-20, -10, 10, 20)
#generate the chains
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k){
  X[i, ] <- cauchy.chain(n, x0[i])
}
  
#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)){
  psi[i,] <- psi[i,] / (1:ncol(psi))
}
  
print(Gelman.Rubin(psi))
#plot psi for the four chains

for (i in 1:k){
  plot(psi[i, (b+1):n], type="l",
       xlab=i, ylab=bquote(psi))
}
  

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n){
  rhat[j] <- Gelman.Rubin(psi[,1:j])
}

plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)
min(which(rhat[(b+1):n]<=1.2))

## -----------------------------------------------------------------------------
set.seed(1)
n<-min(which(rhat[(b+1):n]<=1.2))+b
x<-cauchy.chain(n,-10)
x<-x[-(1:1000)]
a <- ppoints(100)
QR <- qcauchy(a) #quantiles of Cauchy
Q <- quantile(x, a)
qqplot(QR, Q, main="",
       xlab="Cauchy Quantiles", ylab="Sample Quantiles")
hist(x, breaks="scott", main="", xlab="", freq=FALSE)
lines(QR,dcauchy(QR))

a <- c(.05, seq(.1, .9, .1), .95)
Q <- qcauchy(a)
Qx<-quantile(x,a)
print(round(cbind(Q,Qx),3))

## -----------------------------------------------------------------------------
set.seed(3)

bivariate.chain<-function(N,X1,n=15,a=5,b=5){
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  X[1,]<-X1
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, x2)
    x1 <- X[i, 1]
    X[i, 2] <- rbeta(1, x1+a, n-x1+b)
  }
  return(X)
}
k<-4
n<-30000
b<-1000
x0 <- matrix(c(1,3,5,10,0.2,0.4,0.6,0.8),4,2)

Xone <- matrix(0, nrow=k, ncol=n)
for (i in 1:k){
  Xone[i, ] <- bivariate.chain(n, x0[i,])[,1]
}
psione <- t(apply(Xone, 1, cumsum))
for (i in 1:nrow(psione)){
  psione[i,] <- psione[i,] / (1:ncol(psione))
}
print(Gelman.Rubin(psione))
rhatone <- rep(0, n)
for (j in (b+1):n){
  rhatone[j] <- Gelman.Rubin(psione[,1:j])
}

Xtwo <- matrix(0, nrow=k, ncol=n)
for (i in 1:k){
  Xtwo[i, ] <- bivariate.chain(n, x0[i,])[,2]
}
psitwo <- t(apply(Xtwo, 1, cumsum))
for (i in 1:nrow(psitwo)){
  psitwo[i,] <- psitwo[i,] / (1:ncol(psitwo))
}
print(Gelman.Rubin(psitwo))
rhattwo <- rep(0, n)
for (j in (b+1):n){
  rhattwo[j] <- Gelman.Rubin(psitwo[,1:j])
}



plot(rhatone[(b+1):n], type="l", xlab="", ylab="R",col=1,ylim=c(0,2))
lines(rhattwo[(b+1):n],col=2)
abline(h=1.2, lty=2)




## -----------------------------------------------------------------------------
ak<-function(k,d=2,a=c(1,2)){
  result<-sum(a^2)/factorial(k)*(-sum(a^2)/2)^k*(1/(2*k+1)-1/(2*k+2))*exp(lgamma((d+1)/2)+lgamma(k+3/2)-lgamma(k+d/2+1))
  return(result)
}

## -----------------------------------------------------------------------------
Sn<-function(d=2,a=c(1,2),epsilon=1e-7){
k<-1
sn<-0
while (abs(ak(k,d=d,a=a))>=1e-7) {
  k<-k+1
  sn<-sn+ak(k,d=d,a=a)
}
return(sn)
}

## -----------------------------------------------------------------------------
Sn(d=2,a=c(1,2))


## -----------------------------------------------------------------------------
K<-c(4:25,100)
B <- rep(0,length(K))
i<-1
for(k in K){
  ck<-exp(log(2) + lgamma(k/2) - (1/2)*(log(pi*(k-1))) - lgamma((k-1)/2))
  ckk<-exp(log(2) + lgamma((k+1)/2) - (1/2)*(log(pi*(k))) - lgamma((k)/2))
  f<- function(a){
    return((ck*integrate(function(x)(1 + (x^2)/(k-1))^(-k/2), lower = 0, upper = sqrt((a^2)*(k-1)/(k-a^2)))$value) - (ckk*integrate(function(x)(1 + (x^2)/k)^(-(k+1)/2), lower = 0,upper = sqrt((a^2)*k/(k+1-a^2)))$value))
  }
  B[i] <- uniroot(f, lower = 0.01, upper = 1+sqrt(k)/2)$root
  i<-i+1
}


## -----------------------------------------------------------------------------
K<-c(4:25, 100)
BB<-rep(0,length(K))
i<-1
for(k in K){
  ck<-exp(log(2) + lgamma(k/2) - (1/2)*(log(pi*(k-1))) - lgamma((k-1)/2))
  ckk<-exp(log(2) + lgamma((k+1)/2) - (1/2)*(log(pi*(k))) - lgamma((k)/2))
  g <- function(a){
    return((ck*integrate(function(x)(1 + (x^2)/(k-1))^(-k/2),
                           lower = sqrt((a^2)*(k-1)/(k-a^2)), upper = Inf)$value) - (ckk*integrate(function(x)(1 + (x^2)/k)^(-(k+1)/2), 
                                                                                                      lower = sqrt((a^2)*k/(k+1-a^2)), upper = Inf)$value))
  }
  BB[i] <- uniroot(f, lower = 0.01, upper = 1+sqrt(k)/2)$root
  i<-i+1
}
cbind(K,B,BB)


## -----------------------------------------------------------------------------
Y<-c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)

EM<-function(y,max.itr=1000,eps=1e-6){
  i<-1
  n<-length(y)
  m<-n-sum(y==1)
  lambda1<-1
  lambda2<-mean(y)
  while(abs(lambda1-lambda2)>=eps&i<=max.itr){
    lambda1<-lambda2
    lambda2<-(sum(y[1:m])+(n-m)*(1+lambda2))/n
    i<-i+1
  }
  result<-c(lambda2,i)
  return(result)
}

n<-length(Y)
m<-n-sum(Y==1)
(sum(Y[1:m])+(n-m))/m # likelihood
EM(Y)[1] # E.M.


## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)


## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
rsq <- function(mod) summary(mod)$r.squared
lapply(lapply(formulas,lm,data=mtcars),rsq)

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
lapply(lapply(bootstraps,lm,formula=mpg ~ disp),rsq)

## -----------------------------------------------------------------------------
vapply(mtcars,sd,numeric(1))

## -----------------------------------------------------------------------------
df <- data.frame(x = 1:10, y = letters[1:10],z=2:11)
vapply(df[,which(vapply(df, is.numeric,logical(1))==TRUE)],sd,numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  R_Gibbs<-function(n,a,b){
#    N <- 10000 #length of chain
#    X <- matrix(0, N, 2) #the chain, a bivariate sample
#  
#    X[1, ] <- c(1, 0.5) #initialize
#    for (i in 2:N) {
#      x2 <- X[i-1, 2]
#      X[i, 1] <- rbinom(1,n,x2)
#      x1 <- X[i, 1]
#      X[i, 2] <- rbeta(1, x1+a, n-x1+b)
#    }
#  
#    return(X)
#  }
#  
#  
#  library(Rcpp)
#  cppFunction('NumericMatrix C_Gibbs(int n,double a, double b){
#              int N=10000;
#              int x;
#              double y;
#              NumericMatrix M(N,2);
#              M(0,0)=1;
#              M(0,1)=0.5;
#              for(int i=1;i<N;i++){
#              y=M(i-1,1);
#              M(i,0)=rbinom(1,n,y)[0];
#              x=M(i,0);
#              M(i,1)=rbeta(1,x+a,n-x+b)[0];
#              }
#              return(M) ;
#              }')
#  
#  
#  set.seed(0)
#  burn<-1000
#  N<-10000
#  GibbsR=R_Gibbs(15,3,3)
#  GibbsR<-GibbsR[(burn+1):N,]
#  GibbsC=C_Gibbs(15,3,3)
#  GibbsC<-GibbsC[(burn+1):N,]
#  GibbsR_x<-GibbsR[,1]
#  GibbsR_y<-GibbsR[,2]
#  
#  GibbsC_x<-GibbsC[,1]
#  GibbsC_y<-GibbsC[,2]
#  
#  qqplot(GibbsR_x,GibbsC_x)
#  qqplot(GibbsR_y,GibbsC_y)
#  

## ----eval=FALSE---------------------------------------------------------------
#  
#  library(microbenchmark)
#  timee<-microbenchmark(GibbsR=R_Gibbs(15,3,3),GibbsC=C_Gibbs(15,3,3))
#  summary(timee)[, c(1,3,5,6)]
#  

