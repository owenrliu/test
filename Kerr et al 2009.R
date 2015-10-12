#####################################################################
#Reproduction of age-structured model from Kerr et al. 2009
#The role of spatial dynamics in the stability, resilience, and
#productivity of an estuarine fish population
#Owen Liu, owenrliu@gmail.com
#October 2015
#####################################################################

#list of model parameters
n.ages <- 13  # number of age groups
B1 <- 1985000 # maximum number of recruits produced
B2 <- 1000    # controls rate at which the max recruits per spawner is reached
Z1 <- 1.54    # total mortality over the larval period
Za <- 0.59    # mean total adult mortality
Linf <- 217   # asymptotic VB length
kr <- 0.39    # rate (yr^-1) at which the growth model approaches the asymptote
kd <- 0.69    # rate (yr^-1) at which the growth model approaches the asymptote
a0 <- 0       # theoretical age at zero length (VB)
alpha <- 6.42e-6  #length-mass parameter
beta <- 3.28  #length-mass parameter
D <-0.25      #dispersive fraction

####Calculate population of fish in each contingent##################
# 2 different contingents, resident and dispersive (with fraction D)
# separate calculations for age-zero fish (from both contingents
# together, with error), age-1 fish (based on larval mortality, which
# varies between wet and dry years and between contingents), and age
# 2-12 dynamics, which are governed by a universal, constant 'adult'
# mortality rate, Za

N.vec <- function (B1, B2, SSBpop, Z1.vec, Za, D) {
  
  N.vec.r <- numeric(13) #population vector, resident contingent
  
  N.vec.d <- numeric(13) #population vector, dispersive contingent
  
  N0<-(B1*SSBpop)/(B2 + SSBpop) + rnorm(n=1,mean=0,sd=0.03*B1) #age zero
  
  N1r <- N0*(1-D)*exp(-Z1.vec[1]) #Z1.vec[1] is the year-specific larval
                                  #mortality for the resident contigent
  
  N1d <- N0*D*exp(-Z1.vec[2])     #Z1.vec[2] is the year-specific larval
                                  #mortality for the dispersive contigent
  
  N.vec.r[1] <- N0
  N.vec.r[2] <- N1r
  N.vec.d[1] <- N0
  N.vec.d[2] <- N1d
  
  for (i in 2:12) {
    
    N.vec.r[i+1] <- N.vec.r[i]*exp(-Za)
    N.vec.d[i+1] <- N.vec.d[i]*exp(-Za)
    
  }
  
  return(rbind(N.vec.r,N.vec.d))
  
}
test<-N.vec(B1=B1,B2=B2,SSBpop=1000000,Z1.vec=c(1.54,0.5),Za=Za, D=D)
plot(test[1,],type='l')
lines(test[2,],type='l')


####Calculate mass at age for a vector of lengths####
# including VB length-at-age

maa <- function(age.vec,Linf,a0,alpha, beta,k) {
  L.vec<-Linf*(1-exp(-k*(age.vec-a0)))
  M.vec<- alpha*L.vec^beta
  return(M.vec)
}

### Calculate SSB for each contingent####
# as a function of numbers at age at a given time, average
# spawning mass at age, and proportion female spawners at age

calc.SSB <- function (N.vec,M.vec,p.vec) {
  
  SSB.cont <- sum(N.vec*M.vec*p.vec)
  
  return (SSB.cont)
  
}

### Calculate total SSB of the population ####
# calls calc.SSB

SSB.pop <- function (N.vec, M.vec, p.vec) {
  
  calc.SSB (N.vec, M.vec, p.vec)
  
}

####Test simulation####
D <- rnorm(n=100,mean=0.5,sd=0.03) #stochastic dispersive fraction
wet.or.dry<-sample(c(0,1),100,replace=T,prob=c(0.50,0.50)) #wet is zero, dry is 1
drought <- sample(c(0,1),100,replace=T,prob=c(0.20,0.80)) #drought years are zeros
D<-D*drought #dispersive fraction is zero in drought years
Z1r <- wet.or.dry*1.54 #larval mortality, resident contingent
Z1r<-replace(Z1r,Z1r==0,1.2)
Z1d <- wet.or.dry*1.54 #larval mortality, dispersive contingent (env. variable)
Z1d<-replace(Z1d,Z1d==0,0.5)
Z1.vec<-cbind(Z1r,Z1d)

#length, mass, proportion mature
M.vec.r<-maa(age.vec=0:12,Linf=Linf,a0=a0,alpha=alpha, beta=beta,k=kr)
M.vec.d<-maa(age.vec=0:12,Linf=Linf,a0=a0,alpha=alpha, beta=beta,k=kd)
p.vec <- c(0,0,0.73,0.98,1,1,1,1,1,1,1,1,1)

pop.timeseries <- list()
SSB.r <- numeric(100)
SSB.d <- numeric(100)
SSB.total <- numeric(100)

pop.timeseries[[1]] <- N.vec(B1=B1, B2=B2, SSBpop=1E8, Z1.vec=Z1.vec[1,], Za=Za, D=D[1])
SSB.r[1] <- calc.SSB(N.vec=pop.timeseries[[1]][1,],M.vec=M.vec.r,p.vec=p.vec)
SSB.d[1] <- calc.SSB(N.vec=pop.timeseries[[1]][2,],M.vec=M.vec.d,p.vec=p.vec)
SSB.total[1]<-sum(SSB.r[1],SSB.d[1])

for (i in 2:100) {
  pop.timeseries[[i]] <- N.vec(B1=B1, B2=B2, SSBpop=SSB.total[i-1], Z1.vec=Z1.vec[i,], Za=Za, D=D[i])
  SSB.r[i] <- calc.SSB(N.vec=pop.timeseries[[i]][1,],M.vec=M.vec.r,p.vec=p.vec)
  SSB.d[i] <- calc.SSB(N.vec=pop.timeseries[[i]][2,],M.vec=M.vec.d,p.vec=p.vec)
  SSB.total[i]<-sum(SSB.r[i],SSB.d[i])
  
}

#####plots####
plot(1:100,SSB.r/100000,type='l')
plot(1:100,SSB.d/100000,type='l')
plot(1:100,SSB.total/100000,type='l')
plot(pop.timeseries[[1]][1,],type='l')
