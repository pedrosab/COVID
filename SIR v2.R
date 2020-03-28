#SIR model -- SIR (susceptible-infected-removed)
setwd("~/Desktop")
days<-120
Delta_t<-0.001 #in number of days, length of day in a step
steps<-days/Delta_t
t<-seq(from=0,to=steps,by=1)
var_names<-c("S","I","R","Total")
results<-matrix(data=NA, nrow = length(t), ncol = 4, dimnames = list(t,var_names))
beta<-0.4 #infection rate > 0
gamma<-(1/14) #recovery rate >0

#initial conditions, number of people
N<-210e6 #population
Infected<-1
Susceptible<-(N-Infected)
Recovered<-0

#initial conditions, scaled
So<-Susceptible             #So share of susceptible
Io<-Infected                #Share of infected
Ro<-Recovered               #Share of Ro recovered
Sn<-So
In<-Io
Rn<-Ro

for (n in t) {
  results[n+1,]<-c(Sn,In,Rn,Sn+In+Rn)
  Sn_minus_one<-Sn
  In_minus_one<-In
  Rn_minus_one<-Rn
  Sn<-Sn_minus_one+Delta_t*(-beta)*Sn_minus_one/N*In_minus_one
  In<-In_minus_one+Delta_t*(beta*Sn_minus_one/N*In_minus_one-gamma*In_minus_one)
  Rn<-Rn_minus_one+Delta_t*gamma*In_minus_one
}

####################
Max.infected.pop<-max(results[,2])
print(paste("Max share of pop. infected=",round(Max.infected.pop/N,2)))
max.t<-which(results[,2]==Max.infected.pop)-1
print(paste("Day of Max=",round(max.t*Delta_t,2)))

#########################
# estimate of growth rate
Infected.results<-results[,2]
reg<-lm(log(Infected.results[1:max.t])~t[1:max.t])
plot(Infected.results[1:max.t], type="l", main=paste("Growth section of infected curve and forecast with r=",reg$coefficients[2]),col="green")
lines(proj<-exp(reg$coefficients[1]+reg$coefficients[2]*t[1:max.t]), col="black", lty="dashed")
############

print(reg$coefficients[2]/Delta_t)

############
#plot SIR model results
plot(t,results[,1]/N,col="blue", type = "l",main=paste("SIR model with beta=",beta,", gamma=",round(gamma,2),"."), sub="Susceptible=blue, Infected=green, Recovered=red", xlab=paste("Step=",Delta_t,"days"), ylab="Share of Population", ylim = c(0,1))
lines(results[,2]/N,col="green")
lines(results[,3]/N,col="red")
lines(proj/N,col="black", lty="dashed")
############
print(paste("Day 2/Day1-1 % Change",round(results[3/Delta_t,2]/results[2/Delta_t,2]-1,2)))
print(paste("Day 11/Day10-1 % Change",round(results[12/Delta_t,2]/results[11/Delta_t,2]-1,2)))
print(paste("Day 21/Day20-1 % Change",round(results[22/Delta_t,2]/results[21/Delta_t,2]-1,2)))
print(paste("Day 31/Day30-1 % Change",round(results[32/Delta_t,2]/results[31/Delta_t,2]-1,2)))