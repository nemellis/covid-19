library(epimdr)

# SEIR model with vital dynamics for SARS-CoV-2
# assumptions: sterilizing immunity, stable population
seir_ode<-function(t,Y,par){
  S<-Y[1]
  E<-Y[2]
  I<-Y[3]
  R<-Y[4]
  
  beta<-par[1]
  nu<-par[2]
  gamma<-par[3]
  mu<-par[4]
  
  dYdt<-vector(length=4)
  dYdt[1]=mu-beta*I*S-mu*S
  dYdt[2]=beta*I*S-(nu+mu)*E
  dYdt[3]=nu*E-(gamma+mu)*I
  dYdt[4]=gamma*E-mu*R
  
  return(list(dYdt))
}
  
beta<-4.682; # from trajectory mapping
nu<-1/5.1; # Study shows that median incubation is 5.1 days; however, latent period may be shorter 
gamma<-1/14; # WHO says it takes 14 days on avg to recover
mu<-16300/(390144*365) # Orleans Parish birth rate from Nola.com, and population from time-series covid-19 datasheet
init<-c(0.8,0.1,0.1,0) # after National Championship, leading up to Mardi Gras
t<-seq(0,365)
Y<-390144 # from time-series covid-19 datasheet 
par<-c(beta,nu,gamma,mu)
# Solve system using lsoda
sol<-lsoda(init,t,seir_ode,par)

plot(t,sol[,2],type="l",col="blue", xlim=c(0,100), ylim=c(0,1),ylab="Proportion")
lines(t,sol[,3],col="orange")
lines(t,sol[,4],col="red")  
lines(t,1-rowSums(sol[,2:4]),col="green")
legend(60,0.7,legend=c("S","E","I","R"),col=c("blue","orange","red","green"), lty=1, cex=0.8)