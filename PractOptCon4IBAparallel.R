## Rcode to accompany the paper
## A Mathematical Model for the Control of Swimmer's Itch
## J.P. Peirce^{1,3}, J.J. Pellett^1, G.J.  Sandland^2,3
## Code written by J. Peirce (jpeirce@uwlax.edu)
#####
library(deSolve)
library(foreach)
library(doParallel)

###### Practical Optimal Control
m.par <- c(
  rho = 0.745,  # prevalence of parasite in migrating merganser
  s = 0.40, # survival probability of egg to adult merganser
  beta = 1.00667*10^(-8), # transmission rate from infected snail to susceptible merganser
  muB = 1/(10*364.25), # natural death rate of merganser
  kB = 8.632638*10^(-6), # mortality rate in merganser due to infection
  chi = 3.523428*10^(-5), # transmission rate from infected merganser to susceptible snail
  muS = (1-.333)/77, # natural death rate of snail
  gammaB = 1/14, # rate out of exposed bird population
  gammaS = 1/35, # rate out of exposed snail population
  kS = 5.48*10^(-3), # snail death rate due to parasite
  Tf = 7*30  # length of season
)

# Bird Migration Function
LambdaIn <- function(t){
  return( 2.38834*exp(-(t-10.5)^2/55.125) )
}
# Bird Birth Function
bB <- function(t){
  return( 13.733*exp(-(t-42)^2/55.125))
}
# Snail Birth Function
bS <- function(t){
  return( 0.01*sqrt(t)*(0.997)^(t^1.5) )
}

# Initial Conditions
m.ic <- c(
  SB = 0,
  EB = 0,
  IB = 0,
  SS = 1000000,
  ES = 0,
  IS = 0
)
SwimItchModel <-function(t, state, parameters){
  with(as.list(c(state,parameters)), {
    dSB <- (1-rho)*LambdaIn(t) + s*bB(t) - beta*SB*IS - muB*SB
    dEB <- beta*SB*IS - (gammaB+muB)*EB
    dIB <- rho*LambdaIn(t) + gammaB*EB - (muB+kB)*IB
    dSS <- bS(t)*SS - chi*SS*IB - muS*SS
    dES <- chi*SS*IB-(gammaS+muS)*ES
    dIS <- gammaS*ES-(muS+kS)*IS
    list(c(dSB,dEB,dIB,dSS,dES,dIS))
  }) # ends with(as.list())
}
#### time step
dt <- 0.1 # timestep for ode solver

## Treat 1 time function
PracticalModel1 <- function(treat1, max.treat){
  #solve until treatment
  soln.beforeT <- as.data.frame(
    ode(y=m.ic, times=seq(0,treat1, by=dt), func=SwimItchModel, parms=m.par)
  )
  # numbers of at treatment day
  m.ic <- tail(soln.beforeT,1)[-1] # -1 removes the time column
  # new initial conditions
  if (m.ic$SB+m.ic$EB+m.ic$IB < max.treat){ # on that day, all birds are treated
    m.ic$SB <- m.ic$SB + m.ic$EB+m.ic$IB
    m.ic$EB <- 0
    m.ic$IB <- 0
  } else { # if there are more birds present than max.treat, treat % of max.treat
    m.ic$SB <- m.ic$SB + (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat +
      (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$EB <- m.ic$EB - (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$IB <- m.ic$IB - (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
  }
  m.ic <- unlist(m.ic, row)
  #solve after treatment
  soln.afterT <- as.data.frame(
    ode(y=m.ic, times=seq(treat1, m.par["Tf"], by=dt), 
        func=SwimItchModel, parms=m.par)
  )
  return(rbind(soln.beforeT[soln.beforeT$time<treat1,], soln.afterT))
}

start_time <- Sys.time()

max.M <- 80
opt.time1 <- data.frame(
  M = 1:max.M,
  day = numeric(max.M)
)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload computer
registerDoParallel(cl)

# in parallel
for (M in 1:max.M){
  meanIS <- foreach(t=1:(m.par["Tf"]-2), .combine='c') %dopar% {
    library(deSolve)
    m.ic <- c(SB = 0, EB = 0, IB = 0,
              SS = 1000000, ES = 0, IS = 0)
    pract.sol <- PracticalModel1(treat1 = t, max.treat = M)
    # mean only in recreational period
    mean(pract.sol[pract.sol$time >= 60 & pract.sol$time <= 165,]$IS)
  }
  opt.time1[M,"day"] <-which.min(meanIS)
}

## Treat 2 times
# Initial Conditions
m.ic <- c(SB = 0, EB = 0, IB = 0,
          SS = 1000000, ES = 0, IS = 0)

PracticalModel2 <- function(treat1, treat2, max.treat){
  #solve until treatment
  soln.beforeT1 <- as.data.frame(
    ode(y=m.ic, times=seq(0,treat1, by=dt), func=SwimItchModel, parms=m.par)
  )
  # numbers at treatment day
  m.ic <- tail(soln.beforeT1,1)[-1] # -1 removes the time column
  # new initial conditions
  if (m.ic$SB+m.ic$EB+m.ic$IB < max.treat){ # on that day, all birds are treated
    m.ic$SB <- m.ic$SB + m.ic$EB+m.ic$IB
    m.ic$EB <- 0
    m.ic$IB <- 0
  } else { # if there are more birds present than max.treat, treat % of max.treat
    m.ic$SB <- m.ic$SB + (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat +
      (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$EB <- m.ic$EB - (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$IB <- m.ic$IB - (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
  }
  m.ic <- unlist(m.ic, row)
  #solve after treatment
  soln.afterT1 <- as.data.frame(
    ode(y=m.ic, times=seq(treat1, treat2, by=dt), 
        func=SwimItchModel, parms=m.par)
  )
  # new initial conditions
  m.ic <- tail(soln.afterT1,1)[-1] # -1 removes the time column
  # new initial conditions
  if (m.ic$SB+m.ic$EB+m.ic$IB < max.treat){ # on that day, all birds are treated
    m.ic$SB <- m.ic$SB + m.ic$EB+m.ic$IB
    m.ic$EB <- 0
    m.ic$IB <- 0
  } else { # if there are more birds present than max.treat, treat % of max.treat
    m.ic$SB <- m.ic$SB + (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat +
      (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$EB <- m.ic$EB - (m.ic$EB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
    m.ic$IB <- m.ic$IB - (m.ic$IB/(m.ic$SB+m.ic$EB+m.ic$IB))*max.treat
  }
  m.ic <- unlist(m.ic, row)
  #solve after treatment
  soln.afterT2 <- as.data.frame(
    ode(y=m.ic, times=seq(treat2, m.par["Tf"], by=dt), 
        func=SwimItchModel, parms=m.par)
  )
  return(rbind(soln.beforeT1[soln.beforeT1$time < treat1,], 
               soln.afterT1[soln.afterT1$time < treat2,], 
               soln.afterT2))
}

max.M <- 0.5*max.M
opt.time2 <- data.frame(
  M = 1:max.M,
  day1 = numeric(max.M),
  day2 = numeric(max.M)
)

# in parallel
# initial NA matrix for meanIS
meanIS.matrix <-  matrix(data=NA,nrow=m.par["Tf"],ncol=m.par["Tf"])
for (M in 1:max.M){
  for (t1 in 1:(m.par["Tf"]-1)){
    meanIS.matrix[t1,] <- foreach(t2=1:as.numeric(m.par["Tf"]), .combine='c') %dopar% {
      if (t2 > t1 & t2 <= (m.par["Tf"]-1)){
        library(deSolve)
        m.ic <- c(SB = 0, EB = 0, IB = 0,
                  SS = 1000000, ES = 0, IS = 0)
        pract.sol <- PracticalModel2(treat1 = t1, treat2 = t2, max.treat = M)
        mean(pract.sol[pract.sol$time >= 60 & pract.sol$time <= 165,]$IS)
      } else {NA}
    }
  }
  # Location of minimum 
  opt.time2[M,2:3] <-which(meanIS.matrix==min(meanIS.matrix, na.rm = TRUE), arr.ind = TRUE)
}

stopCluster(cl) # closes the parallel cluster of cores

end_time <- Sys.time()
(time <- end_time - start_time)

write.csv(opt.time1, file = "opttime1.csv")
write.csv(opt.time2, file = "opttime2.csv")
