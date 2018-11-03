# Orstein–Uhlenbeck Process
alpha <- c(0.1,1,5)
sigma <- c(0.1,0.2,0.5)
b <- c(-5,5)

## random walk
ran_walk <- function(alpha,sigma,b,r0=1,T=500,delta=1/500) {
  dim<-length(alpha)*length(sigma)*length(b)
  results<-data.frame(matrix(0,(T-r0)/delta,dim)); ran_nrml<-results
  for (i in 1:dim) {ran_nrml[,i]<-rnorm(nrow(ran_nrml))}
  init<-data.frame(matrix(0,3,dim))
  for (j in 1:length(alpha)) {
    for (k in 1:length(sigma)) {
      for (l in 1:length(b)) {
        num <- (j-1)*6+(k-1)*2+l
        init[1,num]<-alpha[j]
        init[2,num]<-sigma[k]
        init[3,num]<-b[l]
      }
    }
  }
  ri<-data.frame(matrix(r0,1,dim))
  for (i in 1:nrow(results)) {
    ri<-exp(-init[1,]*delta)*ri[1,1:dim]+init[3,]*(1-exp(-init[1,]*delta))+
      init[2,]*sqrt((1-exp(-2*init[1,]*delta))/2/init[1,])*ran_nrml[i,]
    results[i,]<-ri[1:dim]
  }
  results
}
results<-ran_walk(alpha,sigma,b)

## plot
time<-seq(1,500,1/500)
results2<-cbind(time,results)
library(ggplot)
for (j in 1:length(alpha)) {
  for (k in 1:length(sigma)) {
    for (l in 1:length(b)) {
      num <- (j-1)*6+(k-1)*2+l
      ggplot() + deom_line(data=results2,aes(x=results2[,1],y=results2[,num])) + 
        labs(x="t",y="r(t)",
             title=expression(paste(alpha,"=",alpha[j],", ",sigma,"=",sigma[k],", b=",b[l]))) + 
        theme(plot.title = element_text(hjust = 0.5))
    }
  }
}

## 

# Poisson Process

## Integration
lamda_t <- function(t) {
  lamdat<- sqrt(t) + exp(-t)*sin(2*pi*t)
  lamdat
}

Z <- integrate(lamda_t,0,5)
Z$value

## poisson process
tau_t<-matrix(0,1,1)
for (i in 1:1000) {
  S<-5* runif(rpois(1,5))
  T<-matrix(0,length(S),1)
for (j in 1: length(S)) {
  if (length(S)==0) next
  lamda_f <- function(t) {
    lamda_f <- 2/3*t^(3/2)-exp(-t)+1-S[j]
    lamda_f
  }
  T[j,1]<-uniroot(lamda_f,c(0,5))$root
  tau<-matrix(0,length(S),1)
}
  for (k in 1:nrow(T)) {
    if (length(T)==0) next
    u0 <- runif(1)
    if (u0<=(sqrt(T[k,1])+exp(-T[k,1])*sin(2*pi*T[k,1]))/(sqrt(T[k,1])+exp(-T[k,1]))) {
      tau[k,1]<-T[k,1]
    }
    else next
  }
  tau_t<-rbind(tau_t,tau)
}
tau_t<-as.data.frame(tau_t[tau_t>0 & tau_t<5])

true_f <- function(t) {
  true_f<- (sqrt(t)+exp(-t)*sin(2*pi*t))/(Z$value)
  true_f
}
library(ggplot2)
ggplot() + stat_function(aes(col="true"), fun=true_f) +
  geom_density(data=tau_t,aes(tau_t[,1], col='kernel')) + 
  scale_color_manual(values = c("red", "blue")) +
  labs(title = paste("kernel density vs. true density with "),
       x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))




