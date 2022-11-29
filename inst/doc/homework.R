## -----------------------------------------------------------------------------
plot(1:10,1:10,xlab="",ylab="")

## -----------------------------------------------------------------------------
data_now <- carData::Freedman
knitr::kable(tail(data_now),format = "html",caption="Table1: Freeedman Dataset")

## -----------------------------------------------------------------------------
table2 <- knitr::kable(tail(data_now),format = "html",caption="Table2: Freeedman Dataset")

table2 <- kableExtra::row_spec(table2,0, color = "white", background = "#696969" ) 
# Set the fill method of the header

table2 <- kableExtra::column_spec(table2,5,color=ifelse(data_now$crime[105:110]>2500,"red","black"),bold=data_now$crime[105:110]>2500)
# Bold the number of crimes greater than 2500 in red

kableExtra::kable_styling(table2,bootstrap_options = "bordered",full_width = F,font_size=20)
# Use bootstrap_options="bordered" to set borders for the table
# "full_width=F" contols the table's width such that it won't fill the whole page
# We can use font_size to set the font size in table

## ---- results='asis'----------------------------------------------------------
cat("Freedman Summary\n")
pander::pander(summary(data_now))
# We can also use pander to show tables.

## -----------------------------------------------------------------------------
set.seed(15) 
u <- runif(n=5000,min=0,max=1) # generate U
pareto <- 2/sqrt(1-u)          # generate Pareto(2,2)

## -----------------------------------------------------------------------------
hist(pareto,probability = T,xlim=c(2,60),breaks = 100,main="Pareto(2,2)",xlab="x")
x <- seq(2,60,length.out=1000)
lines(x,16/x^3,col="red") # line the density of Pareto(2,2)
legend("bottomright",legend = "real density", col="red",lty=1,bty="n")

## -----------------------------------------------------------------------------
# gen_beta generate samples from Beta(a,b) using the acceptance-rejection method
gen_beta <- function(a,b,n){
  # a,b the parameters of B(a,b)
  # n the number of samples
  count <- 0 
  x <- numeric(n)
  c <- 1/beta(a,b)
  while(count < n){
    u <- runif(1)
    y <- runif(1) # random variable from g(.)
    if(y^a*(1-y)^b >= u){
      # accept y
      count<-count+1
      x[count] <- y
    }
  }
  return(x)
}

## -----------------------------------------------------------------------------
set.seed(7)
beta_user <- gen_beta(3,2,1000)

## -----------------------------------------------------------------------------
hist(beta_user,main="Beta(3,2)",probability = T,breaks = 12,xlab="x")
x <- seq(0,1,length.out=1000)
lines(x,x^2*(1-x)/beta(3,2),col="red") # line the density of Beta(3,2)
legend("topleft",legend = "real density", col="red",lty=1,bty="n")

## -----------------------------------------------------------------------------
gen_expgamma <- function(r,b,n){
  # r,b the parameter of lambda
  # n the size of samples
  lambda <- rgamma(n,r,b) # generate lambda from Gamma(r,b) 
  y <- sapply(1:n,function(i) rexp(1,lambda[i])) # generate Y
  return(y)
}

# Generate samples
set.seed(15)
expg_user <- gen_expgamma(4,2,1000)

## -----------------------------------------------------------------------------
hist(expg_user,probability  = T,xlab="x",main="Histogram of Exponential-Gamma mixture")

## -----------------------------------------------------------------------------
gen_beta_byg <- function(a,b,n){
  u <- rgamma(n,a,1)
  v <- rgamma(n,b,1)
  x <- u/(u+v)
  return(x)
}

## -----------------------------------------------------------------------------
n <- 1000 # sample size
rep <- 100 # iterations
a <-b<-2 # parameters

set.seed(15)
t1 <- system.time(for (i in 1:rep)
  gen_beta(a,b,n))
set.seed(15)
t2 <- system.time(for (i in 1:rep)
  gen_beta_byg(a,b,n))
set.seed(15)
t3 <- system.time(for (i in 1:rep)
  rbeta(n,a,b))
m22 <- matrix(c(t1,t2,t3),nrow=3,byrow=T)[,1:3]
rownames(m22) <- c("beta_accept_reject","beta_from_gamma","rbeta")
colnames(m22) <- c("user", "system", "elapsed")
knitr::kable(m22,caption = "a=2,b=2")

## -----------------------------------------------------------------------------
a<-3
b<-2
set.seed(15)
t1 <- system.time(for (i in 1:rep)
  gen_beta(a,b,n))
set.seed(15)
t2 <- system.time(for (i in 1:rep)
  gen_beta_byg(a,b,n))
set.seed(15)
t3 <- system.time(for (i in 1:rep)
  rbeta(n,a,b))

m32 <- matrix(c(t1,t2,t3),nrow=3,byrow=T)[,1:3]
rownames(m32) <- c("beta_accept_reject","beta_from_gamma","rbeta")
colnames(m32) <- c("user", "system", "elapsed")
knitr::kable(m32,caption = "a=3,b=2")

## -----------------------------------------------------------------------------
a<-2
b<-1
set.seed(15)
t1 <- system.time(for (i in 1:rep)
  gen_beta(a,b,n))
set.seed(15)
t2 <- system.time(for (i in 1:rep)
  gen_beta_byg(a,b,n))
set.seed(15)
t3 <- system.time(for (i in 1:rep)
  rbeta(n,a,b))

m21 <- matrix(c(t1,t2,t3),nrow=3,byrow=T)[,1:3]
rownames(m21) <- c("beta_accept_reject","beta_from_gamma","rbeta")
colnames(m21) <- c("user", "system", "elapsed")
knitr::kable(m21,caption = "a=2,b=1")

## -----------------------------------------------------------------------------
quick_sort<-function(x){
  # sorting function
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))} # Recursion
}

sample_num <- c(1e4,2e4,4e4,6e4,8e4)

simu_sort <- function(i){
  # simulation function
  set.seed(i)
  time <- sapply(1:5,function(j) system.time(quick_sort(sample(1:sample_num[j])))[1])
  return(time)
}

# Calculate computation time averaged over 100 simulations
times <- sapply(1:100,simu_sort)
an <- rowMeans(times)
names(an)<-as.character(sample_num)

## -----------------------------------------------------------------------------
tn <- sample_num*log(sample_num)
reg <- lm(an~tn,data=data.frame(tn,an)) # linear regression

# scatter plot and regressino line
plot(tn,an,xlab="nlog(n)",ylab="time:a_n",main="Sorting time-sample numbers")
abline(a=reg$coefficients[1],b=reg$coefficients[2],col="red")

## -----------------------------------------------------------------------------
m <- 10000
set.seed(15)
U <- runif(m)
T1 <- exp(U) # simple MC
T12 <- (exp(1-U)+T1)/2 # antithetic

## -----------------------------------------------------------------------------
mean(T1)

## -----------------------------------------------------------------------------
mean(T12)

## -----------------------------------------------------------------------------
(var(T1)-var(T12))/var(T1)

## -----------------------------------------------------------------------------
g <-function(x) x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)

f1 <- function(x) 1/sqrt(2*pi)*exp(-x^2/2)

f2 <- function(x) x^2*exp(-x)/2

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se<- numeric(2)
fg<-matrix(0,nrow=2,ncol=m)

# using f1
set.seed(15)
x <- rnorm(m)
fg[1,] <- g(x)/f1(x)
theta.hat[1] <-mean(fg[1,])
se[1] <- sd(fg[1,])

# using f2
set.seed(15)
x <- rgamma(m,3,1)
fg[2,] <- g(x)/f2(x)
theta.hat[2] <-mean(fg[2,])
se[2] <- sd(fg[2,])

## -----------------------------------------------------------------------------
rbind(theta.hat,se)

## -----------------------------------------------------------------------------
summary(fg[1,])
summary(fg[2,])

## -----------------------------------------------------------------------------
m <- 10000
g <- function(x) exp(-x)/(1+x^2)*(x>0)*(x<1)
f <- function(x){ exp(-x)/(1-exp(-1))*(x>0)*(x<1)}

# f3, inverse transform method
set.seed(15)
u <- runif(m)
x <- -log(1-u*(1-exp(-1)))
fg <- g(x)/f(x)
theta.im <- mean(fg)
se.im <-sd(fg)

## -----------------------------------------------------------------------------
set.seed(15)
k<-5
n<-m/k
theta_s <- var_s <-numeric(k)
for(i in 1:k){
  u <- runif(n,(i-1)/5,i/5)
  x <- -log(1-(1-exp(-1))*u)
  fg <- g(x)/k/f(x)
  theta_s[i]<-mean(fg)
  var_s[i]<-var(fg)
}

## -----------------------------------------------------------------------------
sum(theta_s)

## -----------------------------------------------------------------------------
sqrt(sum(var_s))

## -----------------------------------------------------------------------------
gen_lognorm <- function(n,mu,sigma){
# n the number of samples
y <- rnorm(n,mu,sigma)
return(exp(y))
}

# generate data
n <-10000
mu<-0
sigma <-1
set.seed(15)
x <- gen_lognorm(n,mu,sigma)

## -----------------------------------------------------------------------------
# calculate CI
x_log <- mean(log(x))
s <- sd(log(x))
cl <- x_log-qt(0.975,n-1)/sqrt(n)*s
cu <- x_log+qt(0.975,n-1)/sqrt(n)*s

## -----------------------------------------------------------------------------
print(cl)

## -----------------------------------------------------------------------------
print(cu)

## -----------------------------------------------------------------------------
# clear up memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
count5test <- function(x,y){
x_center <- x-mean(x)
y_center <- y-mean(y)
outx <- sum(x_center>max(y_center))+sum(x_center<min(y_center))
outy <- sum(y_center>max(x_center))+sum(y_center<min(x_center))
# return 1 (reject) or 0 (not reject H_0)
return(as.integer(max(outx,outy)>5))
}

## -----------------------------------------------------------------------------
m <- 10000
# generate the data and do test
simu_norm <- function(i,n,mu,sigma1,sigma2){
# i the seed
# n the number of samples
# sigma1 standard error for X
# sigma2 standard error for y
set.seed(i)
x <- rnorm(n,mu,sigma1)
y <- rnorm(n,mu,sigma2)
# the result of count5test and F test (with confidence level equals 0.055)
test <-c(count5test(x,y),  as.integer(var.test(x,y)$p.value<0.055))
return(test)
}

## -----------------------------------------------------------------------------
mu<-0
sigma1<-1
sigma2<-1.5

# power for small samples
n <- 20
power_s <- rowMeans(sapply(1:m,function(i) simu_norm(i,n,mu,sigma1,sigma2)))
# power for medium samples
n <- 50
power_m <- rowMeans(sapply(1:m,function(i) simu_norm(i,n,mu,sigma1,sigma2)))
# power for large samples
n <- 100
power_l <- rowMeans(sapply(1:m,function(i) simu_norm(i,n,mu,sigma1,sigma2)))
power_norm <- matrix(c(20,power_s,50,power_m,100,power_l),nrow=3,byrow=T)
colnames(power_norm)<-c("n","count5test","F-test")
knitr::kable(power_norm,caption = "Power for samples from norm distribution")

## -----------------------------------------------------------------------------
# generate the data and do test
simu_gamma <- function(i,n,a1,b1,a2,b2){
# i the seed
# n the number of samples
# a1,a2 shape parameters for X,y
# b1,b2 rate parameters for x,y
set.seed(i)
x <- rgamma(n,a1,b1)
y <- rgamma(n,a2,b2)
# the result of count5test and F test (with confidence level equals 0.055)
test <-c(count5test(x,y),  as.integer(var.test(x,y)$p.value<0.055))
return(test)
}

## -----------------------------------------------------------------------------
a1 <-b1<-1
a2<-b2<-0.5

# power for small samples
n <- 20
power_s <- rowMeans(sapply(1:m,function(i) simu_gamma(i,n,a1,b1,a2,b2)))
# power for medium samples
n <- 50
power_m <- rowMeans(sapply(1:m,function(i) simu_gamma(i,n,a1,b1,a2,b2)))
# power for large samples
n <- 100
power_l <- rowMeans(sapply(1:m,function(i) simu_gamma(i,n,a1,b1,a2,b2)))

power_gamma <- matrix(c(20,power_s,50,power_m,100,power_l),nrow=3,byrow=T)
colnames(power_gamma)<-c("n","count5test","F-test")
knitr::kable(power_gamma,caption = "Power for samples from gamma distribution")

## -----------------------------------------------------------------------------
time_fail <- c(3, 5, 7, 18, 43, 85, 91, 98,100, 130, 230, 487)
lambda_mle <-length(time_fail)/sum(time_fail)
cat("The MLE of the hazard rate is ", round(lambda_mle,5),".")

## -----------------------------------------------------------------------------
B <- 1e4 # times for bootstrap
lambda_b <-numeric(B)
set.seed(99)
# bootstrap
for(i in 1:B){
  lambda_b[i]<-1/mean(sample(time_fail,12,replace = T))
}

# estimate bias and se
lambda_bias <- mean(lambda_b)-lambda_mle
se.lambda <- sd(lambda_b)

## -----------------------------------------------------------------------------
round(lambda_bias,5)

## -----------------------------------------------------------------------------
round(se.lambda,5)

## -----------------------------------------------------------------------------
# clear up memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
library(boot)
time_fail <- c(3, 5, 7, 18, 43, 85, 91, 98,100, 130, 230, 487)
set.seed(76)
R <- 1e4
stat.b <- boot(data=time_fail, statistic=function(x,i) mean(x[i]),R=R)
# calculate CI
ci <- boot.ci(stat.b, type=c("norm","basic","perc","bca"))

## -----------------------------------------------------------------------------
CIs <- matrix(c(ci$normal[2:3],ci$basic[4:5],ci$percent[4:5],ci$bca[4:5]),nrow=4,byrow=T)
colnames(CIs)<-c("lower","upper")
rownames(CIs) <- c("normal","basic","percentile","BCa")

## -----------------------------------------------------------------------------
knitr::kable(CIs,digits=3,caption = "95% CI")

## -----------------------------------------------------------------------------
# compare the length of for CIs
order(CIs[,2]-CIs[,1],decreasing = F)

## -----------------------------------------------------------------------------
# clear up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
m <- 1e4
B <- 199
ci.norm <- ci.basic<-ci.percent <-matrix(NA,m,2)
boot_mean <- function(x,i) mean(x[i])
n<-50 # the number of samples to generate at each time
mu <- 0

# use MC to calculate CI
for(i in 1:m){
  set.seed(i)
  x <- rnorm(n,mu,1) # generate normal distribution with mean=mu and sd=1
  mean.b <- boot(x,statistic=boot_mean,R=B)
  
  # obtain three 95% CIs
  CIs<-boot.ci(mean.b,type=c("norm","basic","perc"))
  ci.norm[i,]<-CIs$norm[2:3]
  ci.basic[i,] <- CIs$basic[4:5]
  ci.percent[i,] <- CIs$percent[4:5]
}

rate <- matrix(nrow=3,ncol=3)
colnames(rate) <- c("coverage","miss on left","miss on right")
rownames(rate)<-c("normal","basic","percentile")

# the porportion of times that the confidence intervals miss on the left
rate[,2]<-c(sum(mu<ci.norm[,1]),sum(mu<ci.basic[,1]),sum(mu<ci.percent[,1]))/m

# the porportion of times that the confidence intervals miss on the right
rate[,3]<-c(sum(mu>ci.norm[,2]),sum(mu>ci.basic[,2]),sum(mu>ci.percent[,2]))/m

# emperical coverage rates
rate[,1]<-1-rate[,2]-rate[,3]

## -----------------------------------------------------------------------------
knitr::kable(rate,digits=3)

## -----------------------------------------------------------------------------
# clear up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
library(bootstrap)
# load the data
data("scor")

## -----------------------------------------------------------------------------
n <- 88
# estimate theta
eigen_now <- eigen(cov(scor)*(n-1)/n)$values
theta.hat <- eigen_now[1]/sum(eigen_now)

# use jackknife to obtain bias and se of theta.hat
theta.jack <- numeric(n)
for(i in 1:n){
  eigen_now <- eigen(cov(scor[-i,])*(n-2)/(n-1))$values
  theta.jack[i] <- eigen_now[1]/sum(eigen_now)
}

bias.jack <- (n-1)*(mean(theta.jack)-theta.hat)
se.jack <- sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
theta <- c(theta.hat,bias.jack,se.jack)
names(theta)<-c("original","bias,jack","se.jack")

## -----------------------------------------------------------------------------
# print results
round(theta,3)

## -----------------------------------------------------------------------------
# clear up memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
library(DAAG)
data("ironslag")

## -----------------------------------------------------------------------------
attach(ironslag)
n <- length(magnetic)
l <- n*(n-1)/2 # possible cases for leave-two-out

e1 <- e2 <- e3 <- e4 <- numeric(l)
ind <-1

for(i in 1:(n-1)){
  for(j in (i+1):n){
    # leave- two-out
    inds <- c(i,j)
    y <- magnetic[-inds]
    x <- chemical[-inds]
    
    # linear
    m1 <- lm(y~x)
    yhat1 <- m1$coef[1]+m1$coef[2]*chemical[inds]
    e1[ind] <- mean((magnetic[inds]-yhat1)^2)
    
    # quadratic
    m2 <- lm(y~ x+I(x^2))
    yhat2 <- m2$coef[1]+m2$coef[2]*chemical[inds]+m2$coef[3]*chemical[inds]^2
    e2[ind]<-mean((magnetic[inds]-yhat2)^2)
    
    # exponential
    m3 <- lm(log(y) ~ x)
    logyhat3 <- m3$coef[1] + m3$coef[2] * chemical[inds]
    yhat3 <- exp(logyhat3)
    e3[ind] <- mean((magnetic[inds] - yhat3)^2)
    
    # log-log
    m4 <- lm(log(y) ~ log(x))
    logyhat4 <- m4$coef[1] + m4$coef[2] * log(chemical[inds])
    yhat4 <- exp(logyhat4)
    e4[ind] <- mean((magnetic[inds] - yhat4)^2)
    
    ind = ind+1
  }
}

## -----------------------------------------------------------------------------
c(mean(e1),mean(e2),mean(e3),mean(e4))

## -----------------------------------------------------------------------------
# clear up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
# generate data
set.seed(129)
x <- rnorm(100)
y <- x^3

## -----------------------------------------------------------------------------
library(StatComp22003)
z<-matrix(c(x,y),ncol=2)

spearman <- function(z,ix){
  x <- z[,1] # leave x as it is
  y <- z[ix,2] #permute rows of y
  return(cor(x,y,method = "spearman"))
}

# permutation
permu <- boot(data=z,statistic = spearman,R=999,sim="permutation")
tb <- c(permu$t0,permu$t)

# p-value for permutation test
p <- mean(abs(tb)>=abs(tb[1]))
p

## -----------------------------------------------------------------------------
cor.test(x,y)

## -----------------------------------------------------------------------------
# clear up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
# generate a chain
rw.laplace <- function(x0, N, sigma){
  # x0 initiial value
  # N the length of the chain
  # sigma sd for proposal distribution
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0 # count acceptance cases
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(abs(x[i-1])-abs(y))){
      x[i] <- y
      k <- k+1
    }else {
      x[i] <- x[i-1]
      
    } 
  }
   return(list(x=x, k=k))
}

# calculate $\hat{R}$ to monitor convergence 
cal.GR <- function(psi){
  # psi[i,j] is the diagnostic statistic psi(X[i,1:j])
  # for chain in i-th row of X
  
  n <- ncol(psi) 
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)  # between variance est.
  psi.w <- apply(psi, 1, var)
  W <- mean(psi.w)        # within variance est.
  v.hat <- W*(n-1)/n+B/n  # var(psi) est.
  r.hat <- sqrt(v.hat/W)  # G-R statistic
  return(r.hat)
 }

## -----------------------------------------------------------------------------
x0 <- c(-10,-5,0,5,10)
sigmas <- c(0.5,1,2,4)
k <- 5
N <- 5000 # length of chains 
b <- 500 # burn-in lenght
X <- list("s1"=matrix(nrow=k,ncol=N),"s2"=matrix(nrow=k,ncol=N),"s3"=matrix(nrow=k,ncol=N),"s4"=matrix(nrow=k,ncol=N)) # collect chains
accept <- matrix(nrow=length(sigmas), ncol=k) # record acceptance rate for each chain


set.seed(225)
for(i in 1:k){
  
  mc <- rw.laplace(x0[i],N,sigmas[1])
  X$s1[i,] <- mc$x
  accept[1,i] <- mc$k/N
  
  mc <- rw.laplace(x0[i],N,sigmas[2])
  X$s2[i,] <- mc$x
  accept[2,i] <- mc$k/N
  
  mc <- rw.laplace(x0[i],N,sigmas[3])
  X$s3[i,] <- mc$x
  accept[3,i] <- mc$k/N
  
  mc <- rw.laplace(x0[i],N,sigmas[4])
  X$s4[i,] <- mc$x
  accept[4,i] <- mc$k/N
}

## -----------------------------------------------------------------------------
psi <- list("s1"=matrix(nrow=k,ncol=N),"s2"=matrix(nrow=k,ncol=N),"s3"=matrix(nrow=k,ncol=N),"s4"=matrix(nrow=k,ncol=N))
R.hat <-list("s1"=numeric(N),"s2"=numeric(N),"s3"=numeric(N),"s4"=numeric(N))
psi$s1 <- t(apply(X$s1,1,cumsum))
psi$s2 <- t(apply(X$s2,1,cumsum))
psi$s3 <- t(apply(X$s3,1,cumsum))
psi$s4 <- t(apply(X$s4,1,cumsum))

for(i in 1: k){
  
  # calculate the statistic psi
  psi$s1[i,] <- psi$s1[i,]/(1:N)
  psi$s2[i,] <- psi$s2[i,]/(1:N)
  psi$s3[i,] <- psi$s3[i,]/(1:N)
  psi$s4[i,] <- psi$s4[i,]/(1:N)
  
  for(j in (b+1):N){
    R.hat$s1[j] <- cal.GR(psi$s1[,1:j])
    R.hat$s2[j] <- cal.GR(psi$s2[,1:j])
    R.hat$s3[j] <- cal.GR(psi$s3[,1:j])
    R.hat$s4[j] <- cal.GR(psi$s4[,1:j])
  }
  
}

## -----------------------------------------------------------------------------
accept <- cbind(sigmas,accept)
colnames(accept) <- c("sigma","x0=-10","x0=-5","x0=0","x0=5","x0=10")
knitr::kable(round(accept,3),caption="Acceptance rates")

## -----------------------------------------------------------------------------
plot(R.hat$s1[(b+1):N],type="l",xlab="",ylab="R",main="sigma=0.5")
abline(h=1.2,col="red")
plot(R.hat$s2[(b+1):N],type="l",xlab="",ylab="R",main="sigma=1")
abline(h=1.2,col="red")
plot(R.hat$s3[(b+1):N],type="l",xlab="",ylab="R",main="sigma=2")
abline(h=1.2,col="red")
plot(R.hat$s4[(b+1):N],type="l",xlab="",ylab="R",main="sigma=4")
abline(h=1.2,col="red")

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
N <- 5000
k <- 4
X <- matrix(nrow=k, ncol=N)
Y <- matrix(nrow=k, ncol=N)
sd <- sqrt(0.19)
rho <- 0.9
Gibbs.gen <- function(x0,y0,N){
  x <- numeric(N)
  y <- numeric(N)
  x[1] <- x0
  y[1] <- y0
  
  u <- runif(N)
  k <- 0 # count acceptance cases
  for (i in 2:N) {
    x[i] <- rnorm(1, rho*y[i-1], sd)
    y[i] <- rnorm(1,rho*x[i],sd)
    
  }
  return(cbind(x,y))
}

# calculate $\hat{R}$ to monitor convergence 
cal.GR <- function(psi){
  # psi[i,j] is the diagnostic statistic psi(X[i,1:j])
  # for chain in i-th row of X
  
  n <- ncol(psi) 
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)
  B <- n*var(psi.means)  # between variance est.
  psi.w <- apply(psi, 1, var)
  W <- mean(psi.w)        # within variance est.
  v.hat <- W*(n-1)/n+B/n  # var(psi) est.
  r.hat <- sqrt(v.hat/W)  # G-R statistic
  return(r.hat)
 }

## -----------------------------------------------------------------------------
x0 <- -2:1
y0 <- 2:-1
# take x0,y0 as initial value, we generate k chains with length N
b <-1000 # burn-in length
psi <- matrix(0,nrow=k,ncol=N)

set.seed(17)
# generate chains
for(i in 1:k){
  xy <- Gibbs.gen(x0[i],y0[i],N)
  X[i,] <- xy[,1]
  Y[i,] <- xy[,2]
  # calculate psi
  psi[i,] <- sapply(1:N, function(j) mean(X[i,1:j]*Y[i,1:j])-mean(X[i,1:j])*mean(Y[i,1:j]))
}



R.hat <- numeric(N)
for(k in (b+1):N){
  R.hat[k] <- cal.GR(psi[,1:k])
}

## -----------------------------------------------------------------------------
plot(x=X[1,-(1:b)],y=Y[1,-(1:b)],xlab="X",ylab="Y")

## -----------------------------------------------------------------------------
# linear regression
x=X[1,-(1:b)]
y=Y[1,-(1:b)]
lmfit <- lm(y~x)

qqnorm(lmfit$residuals)
qqline(lmfit$residuals,col="red")

## -----------------------------------------------------------------------------
var(lmfit$residuals)

## -----------------------------------------------------------------------------
plot(R.hat[-(1:b)],xlab="",ylab="R",type="l")
abline(h=1.2,col="red")

## -----------------------------------------------------------------------------
# clear up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
# function to generate data
gen_xmy <- function(n,alpha,beta,gamma=1){
  # n number of samples
  x <- rnorm(n,5,1)
  em <- rnorm(n,0,1)

  am <- 1
  ay <- 1
  m <- am +alpha*x+em
 
  y <- ay + beta*m + gamma*x + rnorm(n)
  
  return(data.frame("X"=x,"M"=m,"Y"=y))
}

# calculate alpha_hat and beta_hat
cal_a <- function(X,M){
  data <- data.frame(X,M)
  # data is generated by gen_xmy with three components: X,M,Y
  a <- lm(M~X,data=data)$coefficients["X"] # alpha_hat
  
  return(a)
}

cal_b <- function(X,M,Y){
  data <- data.frame(X,M,Y)
  # data is generated by gen_xmy with three components: X,M,Y
  
  b <- lm(Y~X+M,data=data)$coefficients["M"] #beta_hat
  return(b)
}

## -----------------------------------------------------------------------------
n <- c(20,50,100,200)
R<-199

## -----------------------------------------------------------------------------
# We use perm to do permutation test 
perm <- function(x,m,y,R=199){
  ns <- length(x)
  inds <- t(sapply(1:R,function(i) sample(1:ns,ns))) # permute indices
  t0 <- abs(cal_a(x,m)*cal_b(x,m,y))
  tp <- sapply(1:R,function(i) abs(cal_a(x[inds[i,]],m)*cal_b(x,m[inds[i,]],y))) # calculate T' for permutated data
  ts <- c(t0,tp)
  p.value <- mean(ts>t0)
  return(list("pvalue"=p.value,"t0"=t0,"tp"=tp))
}



# simulation fcn 
simu <- function(i,n,alpha,beta,R=199){
  set.seed(i)
  data <- gen_xmy(n,alpha,beta,1)
  p <- perm(data$X,data$M,data$Y,R=R)$pvalue
  return(p)
}

## -----------------------------------------------------------------------------
test_result1 <- matrix(nrow=length(n),ncol=100)
test_result2 <- matrix(nrow=length(n),ncol=100)
test_result3 <- matrix(nrow=length(n),ncol=100)

for(k in 1:length(n)){
  test_result1[k,] <- sapply(1:100,function(i) simu(i,n[k],0,0))<0.05
  test_result2[k,] <- sapply(1:100,function(i) simu(i,n[k],0,1))<0.05
  test_result3[k,] <- sapply(1:100,function(i) simu(i,n[k],1,0))<0.05
}

## -----------------------------------------------------------------------------
m <- matrix(nrow=3,ncol=4)
m[1,]<-rowMeans(test_result1)
m[2,]<-rowMeans(test_result2)
m[3,]<-rowMeans(test_result3)
rownames(m) <- c("scenario 1","scenario 2","scenario 3")
colnames(m) <- as.character(n)
knitr::kable(round(m,5),caption="Type I error")

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
g <- function(alpha,N,b1,b2,b3,f0){
  tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
  p <- 1/(1+tmp)
  mean(p)-f0
}

cal_alpha <- function(N=1e6,b1=0,b2=1,b3=-1,f0){
  uniroot(g,c(-12,0),N=N,b1=b1,b2=b2,b3=b3,f0=f0)$root
}

## -----------------------------------------------------------------------------
set.seed(225)
N<-1e6
x1 <- rpois(N,1)
x2 <- rexp(N,1)
x3 <- rbinom(N,1,0.5)

## -----------------------------------------------------------------------------
f0 <- c(0.1,0.01,0.001,0.0001)
alpha <- sapply(1:4, function(i) cal_alpha(1e6,0,1,-1,f0[i]))

## -----------------------------------------------------------------------------
plot(x=log(f0),y=alpha)

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
dll<-function(lambda,u,v){
  sum((u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v)))
}

## -----------------------------------------------------------------------------
# u[i]<X_i<v[i]
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- c(12,9,28,14,17,1,24,11,25,3)
dll_root <- uniroot(dll,interval = c(0.05,0.5),u=u,v=v,tol=1e-7)

## -----------------------------------------------------------------------------
dll_root$root

## -----------------------------------------------------------------------------
n <- length(u)
lambda <- 0.05 # initial value
lambda1 <- n/(n/lambda+dll(lambda,u,v))
while(abs(lambda1-lambda)>1e-7){
  lambda <- lambda1
  lambda1 <- n/(n/lambda+dll(lambda,u,v))
}

## -----------------------------------------------------------------------------
lambda1

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
a<-list(1,2,3)
typeof(a)

## -----------------------------------------------------------------------------
unlist(a)

## -----------------------------------------------------------------------------
as.vector(a)

## -----------------------------------------------------------------------------
as.vector(a,mode="integer")

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
data.frame(a=integer(),b=logical())

## -----------------------------------------------------------------------------
data.frame(row.names = c("a","b","c"))

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
df <- data.frame(a=c(1,1.5,6,-5),b=c(2,3,5,10))
data.frame(lapply(df,scale01))

## -----------------------------------------------------------------------------
df <- data.frame(a=c(1,1.5,6,-5),b=c(2,3,5,10), c=c("a", "b", "c", "d"))
data.frame(lapply(df,function(x) if(is.numeric(x)) scale01(x) else x))

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
df <- data.frame(a=1:5,b=3:7,c=-1:3)
vapply(df,sd,0)

## -----------------------------------------------------------------------------
df$d <-c("a","b","c","d","e") # now df is a mixed data frame
vapply(df[vapply(df,is.numeric,logical(1))],sd,numeric(1))

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

## -----------------------------------------------------------------------------
library(Rcpp)
gibs <- function(s0,N=1000, b=100, cor=0.9){
  # N the length of chain
  # b burn length
  # s0 initial state
  ss <- matrix(nrow=N,ncol=2)
  ss[1,] <- s0
  sdc <- sqrt(1-cor^2)
  for(i in 2:N){
    ss[i,1] <- rnorm(1,0.9*ss[i-1,2],sdc)
    ss[i,2] <- rnorm(1,0.9*ss[i,1],sdc)
  }
  return(ss)
}

cppFunction('NumericMatrix gib_c(NumericVector s0, int N, int b, double cor) {
  NumericMatrix ss(N,2);
  ss(0,0) = s0[0];
  ss(0,1) = s0[1];
  double sdc = sqrt(1-cor*cor);
  for(int i=1; i<N; i++){
    ss(i,0) = rnorm(1,0.9*ss(i-1,1),sdc)[0];
    ss(i,1) = rnorm(1,0.9*ss(i,0),sdc)[0];
  }
  return ss;
}')

## -----------------------------------------------------------------------------
set.seed(1596)
r_mc <- gibs(c(0,0))

set.seed(1989)
c_mc <- gib_c(c(0,0),1000,100,0.9)

## -----------------------------------------------------------------------------
b <- 100
N<-1000
qqplot(r_mc[(b+1):N,1],c_mc[(b+1):N,1],main="X_t",xlab="X__t in R",ylab="X_t in C")
qqplot(r_mc[(b+1):N,2],c_mc[(b+1):N,2],main="Y_t",xlab="Y__t in R",ylab="Y_t in C")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(gibbsR=gibs(c(0,0)), gibbsC=gib_c(c(0,0),1000,100,0.9))
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
# clean up the memory
rm(list=ls())
gc()

