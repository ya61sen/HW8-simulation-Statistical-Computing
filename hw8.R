# Orstein–Uhlenbeck Process
alpha_init <- c(0.1,1,5)
sigma_init <- c(0.1,0.2,0.5)
b_init <- c(-5,5)

## random walk
ran_walk <- function(alpha,sigma,b,r0=1,T=500,delta=1/500) {
  dim<-length(alpha)*length(sigma)*length(b)
  init<-data.frame(matrix(0,3,dim))
  result_all<-data.frame(matrix(0,(T-0)/delta+1,dim))
  for (j in 1:length(alpha)) {
    for (k in 1:length(sigma)) {
      for (l in 1:length(b)) {
        results<-numeric((T-0)/delta+1);results[1]<-r0; ran_nrml<-rnorm((T-0)/delta)
        num <- (j-1)*6+(k-1)*2+l
        init[1,num]<-alpha[j]
        init[2,num]<-sigma[k]
        init[3,num]<-b[l]
        part_1<-exp(-init[1,]*delta)
        part_2<-init[3,]*(1-exp(-init[1,]*delta))
        part_3<-init[2,]*sqrt((1-exp(-2*init[1,]*delta))/2/init[1,])
        p1<-as.numeric(part_1[num]);p2<-as.numeric(part_2[num]);p3<-as.numeric(part_3[num])
        for (i in 2:length(results)) {
          results[i]<-p1*results[i-1]+p2+p3*ran_nrml[i-1]
        }
        result_all[,num]<-results
      }
    }
  }
  result_all
}
results<-ran_walk(alpha_init,sigma_init,b_init)

## plot
time<-seq(0,500,1/500)
results2<-cbind(time,results)
library(ggplot2)
for (j in 1:length(alpha_init)) {
  for (k in 1:length(sigma_init)) {
    for (l in 1:length(b_init)) {
      num <- (j-1)*6+(k-1)*2+l
      print(ggplot() + geom_line(data=results2,aes(x=results2[,1],y=results2[,num+1])) + 
        labs(x="t",y="r(t)", title=paste("alpha =",paste(alpha_init[j]),
                      ", sigma =",paste(sigma_init[k]),", b =",paste(b_init[l]))) + 
        theme(plot.title = element_text(hjust = 0.5)))
    }
  }
}

## Euler–Maruyama method
deltas <- c(1,0.5,0.1,0.01)

em_mthd <- function(delta, N, r0=1, alpha=5, sigma=0.1, b=5) {
  sample_N <- NULL
  for (j in 1:N) {
    rt<-numeric(1/delta+1);rt[1]<-r0;r_nrml<-rnorm(1/delta)
    for (i in 2:length(rt)) {
      mu<-rt[i-1]+alpha*(b-rt[i-1])*delta
      sd<-sigma*delta
      rt[i]<-mu+sd*r_nrml[i-1]
    }
    sample_N<-data.frame(rbind(sample_N,rt[length(rt)]))
  }
  sample_N
}

sample_1000 <- data.frame(matrix(0,1000,length(deltas)))
for (i in 1: length(deltas)) {
  sample_1000[,i] <- em_mthd(deltas[i],1000)
}

true_f <- function(x) {
  dnorm(x,(5-4*exp(-5)),sqrt((1-exp(-10))/1000))
}
library(ggplot2)
ggplot()+ stat_function(aes(col="true"),fun=true_f)+
    geom_density(data=sample_1000,aes(sample_1000[,1], col='kernel')) + 
    scale_color_manual(values = c("red", "blue")) + 
    xlim(4,24) + labs(title = paste("delta =",deltas[1]),
         x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))
ggplot()+ stat_function(aes(col="true"),fun=true_f)+
  geom_density(data=sample_1000,aes(sample_1000[,2], col='kernel')) + 
  scale_color_manual(values = c("red", "blue")) + 
  xlim(-5,6) + labs(title = paste("delta =",deltas[2]),
                    x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))
ggplot()+ stat_function(aes(col="true"),fun=true_f)+
  geom_density(data=sample_1000,aes(sample_1000[,3], col='kernel')) + 
  scale_color_manual(values = c("red", "blue")) + 
  xlim(4,6) + labs(title = paste("delta =",deltas[3]),
                    x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))
ggplot()+ stat_function(aes(col="true"),fun=true_f)+
  geom_density(data=sample_1000,aes(sample_1000[,4], col='kernel')) + 
  scale_color_manual(values = c("red", "blue")) + 
  xlim(4,6) + labs(title = paste("delta =",deltas[4]),
                    x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))

# Poisson Process

## Integration
lamda_t <- function(t) {
  lamdat<- sqrt(t) + exp(-t)*sin(2*pi*t)
  lamdat
}
Z <- integrate(lamda_t,0,5)
Z$value

## poisson process
samples <- function(N) {
  int0_5<-2/3*5^(3/2)-exp(-5)+1
  tau_t<-matrix(0,1,1)
  for (i in 1:N) {
    S<-int0_5* runif(rpois(1,int0_5))
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
  tau_t
}

tau_t <- samples(1000)

true_f <- function(t) {
  true_f<- (sqrt(t)+exp(-t)*sin(2*pi*t))/(Z$value)
  true_f
}
library(ggplot2)
ggplot() + stat_function(aes(col="true"), fun=true_f) +
  geom_density(data=tau_t,aes(tau_t[,1], col='kernel')) + 
  scale_color_manual(values = c("red", "blue")) + xlim(0,5) +
  labs(title = paste("kernel density vs. true density"),
       x = "t", y = "value") + theme(plot.title = element_text(hjust = 0.5))

