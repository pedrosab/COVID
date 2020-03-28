#SIR model -- SIR (susceptible-infected-removed)
setwd("~/Desktop")
days<-60
Delta_t<-0.0001 #in number of days, length of day in a step
steps<-days/Delta_t
t<-seq(1:steps)
var_names<-c("S","I","R","Total")
results<-matrix(data=NA, nrow = length(t), ncol = 4, dimnames = list(t,var_names))
beta<-2.2 #infection rate > 0
gamma<-(1/14) #recovery rate >0

#initial conditions, number of people
N<-250000000 #population
Infected<-1
Susceptible<-(N-Infected)
Recovered<-0

#initial conditions, scaled
So<- Susceptible/N            #So share of susceptible
Io<-Infected/N      #Share of infected
Ro<-Recovered/N               #Share of Ro recovered
Sn<-So
In<-Io
Rn<-Ro


for (n in t) {
  Sn_minus_one<-Sn
  In_minus_one<-In
  Rn_minus_one<-Rn
  Sn<-Sn_minus_one+Delta_t*(-beta)*Sn_minus_one*In_minus_one #S(n) as a function of S(n-1)
  In<-In_minus_one+Delta_t*(beta*Sn_minus_one*In_minus_one-gamma*In_minus_one)
  Rn<-Rn_minus_one+Delta_t*gamma*In_minus_one
  results[n,1]<-Sn
  results[n,2]<-In
  results[n,3]<-Rn
  results[n,4]<-Sn+In+Rn
}
plot(t,results[,1],col="blue", type = "l",main=paste("SIR model with beta=",beta,", gamma=",round(gamma,2),"."), sub="Susceptible=blue, Infected=green, Recovered=red", xlab=paste("Step=",Delta_t,"days"), ylab="Share of popultion",ylim=c(0,1))
lines(results[,2],col="green")
lines(results[,3],col="red")
#plot(results[,4])

