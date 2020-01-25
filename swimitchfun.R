## Rcode to accompany the paper
## A Mathematical Model for the Control of Swimmer's Itch
## J.P. Peirce, J.J. Pellett, G.J. Sandland
## Code written by J. Peirce (jpeirce@uwlax.edu)
#####

### Forward-backward sweep method with a 4th order Runge-Kutta scheme to solve the state equations 
# and their corresponding adjoint equations. The forward-backward sweep method takes an initial 
# estimate for the control u and then solves the state equations forward in time using the 
# Runge-Kutta method with the initial conditions. Then, using the state values, the adjoint 
# equations are solved backwards in time using the Runge-Kutta method with the transversality 
# conditions. At this point, the optimal control is updated using the characterization for the 
# optimal control from equation and the values for the state and adjoint variables. This updated 
# control replaces the initial control and the process is repeated until the successive iterates of 
# control values are suffciently close. Motified from method code provided by 
# Lenhart and Workman.

SwimItchControl <- function(m.par, m.ic, A, umax){
  with(as.list(c(m.par,m.ic)), {

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
    
test <- -1 # convergence test variable, 
  # begins the while loop once convergence occurs, 
  # test will become non-negative and while loop ends
    
delta <- 0.01 # tolerance for convergence...when we should say that solution 
  # close enough to optimal
M <- 1000 # number of time steps to use
t <- seq(0,Tf,length=M+1) # creates the time vector from 0 to Tf with M steps
h <- Tf/M # length of one time step
h2 <- h/2 # half of time step used in Runge-Kutta method below
    
# Initializes the population vectors to zero  
SB <- matrix(0,1,M+1)
EB <- matrix(0,1,M+1)
IB <- matrix(0,1,M+1)
SS <- matrix(0,1,M+1)
ES <- matrix(0,1,M+1)
IS <- matrix(0,1,M+1)
    
# Sets first value to the initial population sizes
SB[1] <- m.ic["SB"]
EB[1] <- m.ic["EB"]
IB[1] <- m.ic["IB"]
SS[1] <- m.ic["SS"]
ES[1] <- m.ic["EB"]
IS[1] <- m.ic["IB"]
    
# Initializes the control adjoint vectors to zero  
lambda1 <- matrix(0,1,M+1)
lambda2 <- matrix(0,1,M+1)
lambda3 <- matrix(0,1,M+1)
lambda4 <- matrix(0,1,M+1)
lambda5 <- matrix(0,1,M+1)
lambda6 <- matrix(0,1,M+1)
    
# Initalizes treatment to 0 
u <- matrix(0,1,M+1)
    
# A counter for the number of iterations before convergence
iter <- 0
    
while (test < 0){
  # Used for comparison later to determine if convergence criteria is met
  oldu <- u
  oldSB <- SB
  oldEB <- EB
  oldIB <- IB
  oldSS <- SS
  oldES <- ES
  oldIS <- IS
  oldlambda1 <- lambda1
  oldlambda2 <- lambda2
  oldlambda3 <- lambda3
  oldlambda4 <- lambda4
  oldlambda5 <- lambda5
  oldlambda6 <- lambda6
      
  iter <- iter+1
# Displays count and value of test value....more for my curiosity
# print(iter)
# print(test)
      
# This is an application of a 4th-order Runge-Kutta method used to 
# numericall approximate solutions to differential equations
# See Lenhert and Workman for a good introduction
      
  for (i in 1:M) {
        
    if ((SB[i]+EB[i]+IB[i])==0){
      k11 <- (1-rho)*LambdaIn(t[i])+ s*bB(t[i])-beta*SB[i]*IS[i]-muB*SB[i]
      k12 <- beta*SB[i]*IS[i]-(gammaB+muB)*EB[i]
      k13 <- rho*LambdaIn(t[i])+gammaB*EB[i]-(muB+kB)*IB[i]
      } else {
        k11 <- (1-rho)*LambdaIn(t[i])+ s*bB(t[i])-beta*SB[i]*IS[i]-muB*SB[i]+
          u[i]*(EB[i]+IB[i])/(SB[i]+EB[i]+IB[i])
        k12 <- beta*SB[i]*IS[i]-(gammaB+muB)*EB[i] - 
          u[i]*EB[i]/(SB[i]+EB[i]+IB[i]) 
        k13 <- rho*LambdaIn(t[i])+gammaB*EB[i]-(muB+kB)*IB[i] -
          u[i]*IB[i]/(SB[i]+EB[i]+IB[i])
      }
    k14 <- bS(t[i])*SS[i] - chi*SS[i]*IB[i]-muS*SS[i]
    k15 <- chi*SS[i]*IB[i]-(gammaS+muS)*ES[i]
    k16 <- gammaS*ES[i]-(muS+kS)*IS[i]
    
    if ((SB[i]+h2*k11+EB[i]+h2*k12+IB[i]+h2*k13)==0){
      k21 <- (1-rho)*LambdaIn(t[i]+h2)+ s*bB(t[i]+h2)-beta*(SB[i]+h2*k11)*(IS[i]+h2*k16)-
            muB*(SB[i]+h2*k11)
      k22 <- beta*(SB[i]+h2*k11)*(IS[i]+h2*k16)-(gammaB+muB)*(EB[i]+h2*k12)
      k23 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k12)-(muB+kB)*(IB[i]+h2*k13)      
      } else {
        k21 <- (1-rho)*LambdaIn(t[i]+h2)+ s*bB(t[i]+h2)-beta*(SB[i]+h2*k11)*(IS[i]+h2*k16)-
            muB*(SB[i]+h2*k11) +  
            0.5*(u[i]+u[i+1])*(EB[i]+h2*k12+IB[i]+h2*k13)/(SB[i]+h2*k11+EB[i]+h2*k12+IB[i]+h2*k13)
        k22 <- beta*(SB[i]+h2*k11)*(IS[i]+h2*k16)-(gammaB+muB)*(EB[i]+h2*k12) - 
            0.5*(u[i]+u[i+1])*(EB[i]+h2*k12)/(SB[i]+h2*k11+EB[i]+h2*k12+IB[i]+h2*k13)
        k23 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k12)-(muB+kB)*(IB[i]+h2*k13) -
            0.5*(u[i]+u[i+1])*(IB[i]+h2*k13)/(SB[i]+h2*k11+EB[i]+h2*k12+IB[i]+h2*k13)
        }
    k24 <- bS(t[i]+h2)*(SS[i]+h2*k14) - chi*(SS[i]+h2*k14)*(IB[i]+h2*k13)-muS*(SS[i]+h2*k14)
    k25 <- chi*(SS[i]+h2*k14)*(IB[i]+h2*k13)-(gammaS+muS)*(ES[i]+h2*k15)
    k26 <- gammaS*(ES[i]+h2*k15)-(muS+kS)*(IS[i]+h2*k16)

    if ((SB[i]+h2*k21+EB[i]+h2*k22+IB[i]+h2*k23)==0){
      k31 <- (1-rho)*LambdaIn(t[i]+h2)+s*bB(t[i]+h2)-
        beta*(SB[i]+h2*k21)*(IS[i]+h2*k26)-muB*(SB[i]+h2*k21)
      k32 <- beta*(SB[i]+h2*k21)*(IS[i]+h2*k26)-(gammaB+muB)*(EB[i]+h2*k22)
      k33 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k22)-(muB+kB)*(IB[i]+h2*k23)
      } else {
        k31 <- (1-rho)*LambdaIn(t[i]+h2)+ s*bB(t[i]+h2)-beta*(SB[i]+h2*k21)*(IS[i]+
            h2*k26)-muB*(SB[i]+h2*k21) + 
            0.5*(u[i]+u[i+1])*(EB[i]+h2*k22+IB[i]+h2*k23)/(SB[i]+h2*k21+EB[i]+h2*k22+IB[i]+h2*k23)
        k32 <- beta*(SB[i]+h2*k21)*(IS[i]+h2*k26)-(gammaB+muB)*(EB[i]+h2*k22) - 
            0.5*(u[i]+u[i+1])*(EB[i]+h2*k22)/(SB[i]+h2*k21+EB[i]+h2*k22+IB[i]+h2*k23)
        k33 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k22)-(muB+kB)*(IB[i]+h2*k23) -
            0.5*(u[i]+u[i+1])*(IB[i]+h2*k23)/(SB[i]+h2*k21+EB[i]+h2*k22+IB[i]+h2*k23)
        }
    k34 <- bS(t[i]+h2)*(SS[i]+h2*k24) - chi*(SS[i]+h2*k24)*(IB[i]+h2*k23)-muS*(SS[i]+h2*k24)
    k35 <- chi*(SS[i]+h2*k24)*(IB[i]+h2*k23)-(gammaS+muS)*(ES[i]+h2*k25)
    k36 <- gammaS*(ES[i]+h2*k25)-(muS+kS)*(IS[i]+h2*k26)
        
    if ((SB[i]+h2*k31+EB[i]+h2*k32+IB[i]+h2*k33)==0){
      k41 <- (1-rho)*LambdaIn(t[i+1])+s*bB(t[i+1])-
        beta*(SB[i]+h2*k31)*(IS[i]+h2*k36)-muB*(SB[i]+h2*k31)
      k42 <- beta*(SB[i]+h2*k33)*(IS[i]+h2*k36)-(gammaB+muB)*(EB[i]+h2*k32)
      k43 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k32)-(muB+kB)*(IB[i]+h2*k33)
      } else {
        k41 <- (1-rho)*LambdaIn(t[i+1])+ s*bB(t[i+1])-beta*(SB[i]+h2*k31)*(IS[i]+h2*k36)-
            muB*(SB[i]+h2*k31) + 
            u[i+1]*(EB[i]+h2*k32+IB[i]+h2*k33)/(SB[i]+h2*k31+EB[i]+h2*k32+IB[i]+h2*k33)
        k42 <- beta*(SB[i]+h2*k33)*(IS[i]+h2*k36)-(gammaB+muB)*(EB[i]+h2*k32) -
            u[i+1]*(EB[i]+h2*k32)/(SB[i]+h2*k31+EB[i]+h2*k32+IB[i]+h2*k33)
        k43 <- rho*LambdaIn(t[i]+h2)+gammaB*(EB[i]+h2*k32)-(muB+kB)*(IB[i]+h2*k33) - 
            u[i+1]*(IB[i]+h2*k33)/(SB[i]+h2*k31+EB[i]+h2*k32+IB[i]+h2*k33)
        }
    k44 <- bS(t[i+1])*(SS[i]+h2*k34) - chi*(SS[i]+h2*k34)*(IB[i]+h2*k33)-muS*(SS[i]+h2*k34)
    k45 <- chi*(SS[i]+h2*k34)*(IB[i]+h2*k33)-(gammaS+muS)*(ES[i]+h2*k35)
    k46 <- gammaS*(ES[i]+h2*k35)-(muS+kS)*(IS[i]+h2*k36)
        
    SB[i+1] <- SB[i]+(h/6)*(k11+2*k21+2*k31+k41)
    EB[i+1] <- EB[i]+(h/6)*(k12+2*k22+2*k32+k42)
    IB[i+1] <- IB[i]+(h/6)*(k13+2*k23+2*k33+k43)
    SS[i+1] <- SS[i]+(h/6)*(k14+2*k24+2*k34+k44)
    ES[i+1] <- ES[i]+(h/6)*(k15+2*k25+2*k35+k45)
    IS[i+1] <- IS[i]+(h/6)*(k16+2*k26+2*k36+k46)
  }

for (i in 1:M) {
  j <- M+2-i
 
  if ((SB[j]+EB[j]+IB[j])==0){
    k11 <- (beta*IS[j]+muB)*lambda1[j]-(beta*IS[j])*lambda2[j]
    k12 <- (gammaB+muB)*lambda2[j] - gammaB*lambda3[j]
    k13 <- (muB+kB)*lambda3[j]+chi*SS[j]*lambda4[j]-chi*SS[j]*lambda5[j]
    } else {
    k11 <- (beta*IS[j]+muB+u[j]*(EB[j]+IB[j])/(SB[j]+EB[j]+IB[j])^2)*lambda1[j]-
      (beta*IS[j]+u[j]*EB[j]/(SB[j]+EB[j]+IB[j])^2)*lambda2[j] - 
      u[j]*IB[j]/(SB[j]+EB[j]+IB[j])^2*lambda3[j]
    k12 <- - u[j]*SB[j]/(SB[j]+EB[j]+IB[j])^2*lambda1[j] + (gammaB+muB+
      u[j]*(SB[j]+IB[j])/(SB[j]+EB[j]+IB[j])^2)*lambda2[j] - (gammaB+u[j]*IB[j]/(SB[j]+EB[j]+IB[j])^2)*lambda3[j]
    k13 <- - u[j]*SB[j]/(SB[j]+EB[j]+IB[j])^2*lambda1[j] - u[j]*EB[j]/(SB[j]+EB[j]+IB[j])^2*lambda2[j]+
      (muB+kB+u[j]*(SB[j]+EB[j])/(SB[j]+EB[j]+IB[j])^2)*lambda3[j]+chi*SS[j]*lambda4[j]-
      chi*SS[j]*lambda5[j]
    }
  k14 <- (-bS(t[j])+chi*IB[j]+muS)*lambda4[j]-chi*IB[j]*lambda5[j]
  k15 <- (gammaS+muS)*lambda5[j] - gammaS*lambda6[j]
  k16 <- -A + beta*SB[j]*lambda1[j] - beta*SB[j]*lambda2[j]+(muS+kS)*lambda6[j]
                                     
  if ((0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))== 0){
    k21 <- (beta*0.5*(IS[j]+IS[j-1])+muB)*(lambda1[j]-h2*k11)-
      (beta*0.5*(IS[j]+IS[j-1]))*(lambda2[j]-h2*k12)
    k22 <- (gammaB+muB)*(lambda2[j]-h2*k12) - gammaB*(lambda3[j]-h2*k13)
    k23 <- (muB+kB)*(lambda3[j]-h2*k13)+chi*0.5*(SS[j]+SS[j-1])*(lambda4[j]-h2*k14)-
      chi*0.5*(SS[j]+SS[j-1])*(lambda5[j]-h2*k15)
  } else {
    k21 <- (beta*0.5*(IS[j]+IS[j-1])+muB+0.5*(u[j]+u[j-1])* 
              (0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))/(0.5*(SB[j]+SB[j-1])+
              0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2)*(lambda1[j]-h2*k11)-
      (beta*0.5*(IS[j]+IS[j-1])+0.5*(u[j]+u[j-1])*0.5*(EB[j]+EB[j-1])/(0.5*(SB[j]+SB[j-1])+
              0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2)*(lambda2[j] -h2*k12) -
    0.5*(u[j]+u[j-1])*0.5*(IB[j]+IB[j-1])/(0.5*(SB[j]+SB[j-1])+
              0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda3[j]-h2*k13)
    k22 <- -  0.5*(u[j]+u[j-1])*0.5*(SB[j]+SB[j-1])/(0.5*(SB[j]+SB[j-1])+
              0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda1[j]-h2*k11) + 
              (gammaB+muB+0.5*(u[j]+u[j-1])*(0.5*(SB[j]+SB[j-1]+
              0.5*(IB[j]+IB[j-1])))/(0.5*(SB[j]+SB[j-1])+
              0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2)*(lambda2[j]-h2*k12) 
              (gammaB+0.5*(u[j]+u[j-1])*0.5*(IB[j]+IB[j-1])/(0.5*(SB[j]+SB[j-1])+
             0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2)*(lambda3[j]-h2*k13)
    k23 <- - 0.5*(u[j]+u[j-1])*0.5*(SB[j]+SB[j-1])/(0.5*(SB[j]+SB[j-1])+
             0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda1[j]-h2*k11) - 
      0.5*(u[j]+u[j-1])*0.5*(EB[j]+EB[j-1])/(0.5*(SB[j]+SB[j-1])+
             0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda2[j]-h2*k12)+
      (muB+kB+0.5*(u[j]+u[j-1])*(0.5*(SB[j]+SB[j-1]+
      0.5*(EB[j]+EB[j-1])))/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
          0.5*(IB[j]+IB[j-1]))^2)*(lambda3[j]-h2*k13)+
      chi*0.5*(SB[j]+SB[j-1])*(lambda4[j]-h2*k14)-chi*0.5*(SS[j]+SS[j-1])*(lambda5[j]-h2*k15)
  }
  k24 <- (-bS(t[j-1])+chi*0.5*(IB[j]+IB[j-1])+muS)*(lambda4[j]-h2*k14)-
    chi*0.5*(IB[j]+IB[j-1])*(lambda5[j]-h2*k15)
  k25 <- (gammaS+muS)*(lambda5[j]-h2*k15) - gammaS*(lambda6[j]-h2*k16)
  k26 <- -A + beta*0.5*(SB[j]+SB[j-1])*(lambda1[j]-h2*k11) - 
    beta*0.5*(SB[j]+SB[j-1])*(lambda2[j]-h2*k12)+(muS+kS)*(lambda6[j]-h2*k16)

  if ((0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))==0){
    k31 <- (beta*0.5*(IS[j]+IS[j-1])+muB)*(lambda1[j]-h2*k21)-
      (beta*0.5*(IS[j]+IS[j-1]))*(lambda2[j]-h2*k22)
    k32 <- (gammaB+muB)*(lambda2[j]-h2*k22) - gammaB*(lambda3[j]-h2*k23)
    k33 <- (muB+kB)*(lambda3[j]-h2*k23)+chi*0.5*(SS[j]+SS[j-1])*(lambda4[j]-h2*k24)-
      chi*0.5*(SS[j]+SS[j-1])*(lambda5[j]-h2*k25)
  } else {
    k31 <- (beta*0.5*(IS[j]+IS[j-1])+muB+0.5*(u[j]+u[j-1])*(0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))^2)*(lambda1[j]-h2*k21)-(beta*0.5*(IS[j]+IS[j-1])+
      0.5*(u[j]+u[j-1])*0.5*(EB[j]+EB[j-1])/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))^2)*(lambda2[j] -h2*k22)- 0.5*(u[j]+u[j-1])*
      0.5*(IB[j]+IB[j-1])/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))^2*(lambda3[j]-h2*k23)
    k32 <- -  0.5*(u[j]+u[j-1])*0.5*(SB[j]+SB[j-1])/(0.5*(SB[j]+SB[j-1])+
      0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda1[j]-h2*k21) + 
      (gammaB+muB+0.5*(u[j]+u[j-1])*(0.5*(SB[j]+SB[j-1]+
      0.5*(IB[j]+IB[j-1])))/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))^2)*(lambda2[j]-h2*k22) 
      (gammaB+0.5*(u[j]+u[j-1])*0.5*(IB[j]+IB[j-1])/(0.5*(SB[j]+SB[j-1])+
      0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2)*(lambda3[j]-h2*k23)
    k33 <- - 0.5*(u[j]+u[j-1])*0.5*(SB[j]+SB[j-1])/(0.5*(SB[j]+SB[j-1])+
      0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda1[j]-h2*k21) - 
      0.5*(u[j]+u[j-1])*0.5*(EB[j]+EB[j-1])/(0.5*(SB[j]+SB[j-1])+
      0.5*(EB[j]+EB[j-1])+0.5*(IB[j]+IB[j-1]))^2*(lambda2[j]-h2*k22)+
      (muB+kB+0.5*(u[j]+u[j-1])*(0.5*(SB[j]+SB[j-1]+
      0.5*(EB[j]+EB[j-1])))/(0.5*(SB[j]+SB[j-1])+0.5*(EB[j]+EB[j-1])+
      0.5*(IB[j]+IB[j-1]))^2)*(lambda3[j]-h2*k23)+
      chi*0.5*(SB[j]+SB[j-1])*(lambda4[j]-h2*k24)-chi*0.5*(SS[j]+SS[j-1])*(lambda5[j]-h2*k25)
  }
  k34 <- (-bS(t[j]-h2)+chi*0.5*(IB[j]+IB[j-1])+muS)*(lambda4[j]-h2*k24)-
    chi*0.5*(IB[j]+IB[j-1])*(lambda5[j]-h2*k25)
  k35 <- (gammaS+muS)*(lambda5[j]-h2*k25) - gammaS*(lambda6[j]-h2*k26)
  k36 <- -A + beta*0.5*(SB[j]+SB[j-1])*(lambda1[j]-h2*k21) - 
    beta*0.5*(SB[j]+SB[j-1])*(lambda2[j]-h2*k22)+(muS+kS)*(lambda6[j]-h2*k26)
  
  if ((SB[j-1]+EB[j-1]+IB[j-1])==0){
    k41 <- (beta*IS[j-1]+muB)*(lambda1[j]-h2*k31)-(beta*IS[j-1])*(lambda2[j]-h2*k32)
    k42 <- (gammaB+muB)*(lambda2[j]-h2*k32) - gammaB*(lambda3[j]-h2*k33)
    k43 <- (muB+kB)*(lambda3[j]-h2*k33)+chi*SS[j-1]*(lambda4[j]-h2*k34)-
      chi*SS[j-1]*(lambda5[j]-h2*k35)
  } else {
    k41 <- (beta*IS[j-1]+muB+0.5*u[j-1]*(EB[j-1]+IB[j-1])/(SB[j-1]+EB[j-1]+IB[j-1])^2)*(lambda1[j]-h2*k31)-
      (beta*IS[j-1]+ u[j-1]*EB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2)*(lambda2[j] -h2*k32)-
      u[j-1]*IB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2*(lambda3[j]-h2*k33)
    k42 <- -  u[j-1]*SB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2*(lambda1[j]-h2*k31) + 
      (gammaB+muB+u[j-1]*(SB[j-1]+IB[j-1])/(SB[j-1]+EB[j-1]+IB[j-1])^2)*(lambda2[j]-h2*k32)-
      (gammaB+u[j-1]*IB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2)*(lambda3[j]-h2*k33)
    k43 <- -u[j-1]*SB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2*(lambda1[j]-h2*k31) - 
      u[j-1]*EB[j-1]/(SB[j-1]+EB[j-1]+IB[j-1])^2*(lambda2[j]-h2*k32)+
      (muB+kB+u[j-1]*(SB[j-1]+EB[j-1])/(SB[j-1]+EB[j-1]+IB[j-1])^2)*(lambda3[j]-h2*k33)+
      chi*SB[j-1]*(lambda4[j]-h2*k34)-chi*SS[j-1]*(lambda5[j]-h2*k35)
  }
  k44 <- (-bS(t[j-1])+chi*IB[j-1]+muS)*(lambda4[j]-h2*k34)-
    chi*IB[j-1]*(lambda5[j]-h2*k35)
  k45 <- (gammaS+muS)*(lambda5[j]-h2*k35) - gammaS*(lambda6[j]-h2*k36)
  k46 <- -A + beta*SB[j-1]*(lambda1[j]-h2*k31) - 
    beta*SB[j-1]*(lambda2[j]-h2*k32)+(muS+kS)*(lambda6[j]-h2*k36)
  
  lambda1[j-1]  <-  lambda1[j] - (h/6)*(k11 + 2*k21 + 2*k31 + k41)
  lambda2[j-1]  <-  lambda2[j] - (h/6)*(k12 + 2*k22 + 2*k32 + k42)
  lambda3[j-1]  <-  lambda3[j] - (h/6)*(k13 + 2*k23 + 2*k33 + k43)
  lambda4[j-1]  <-  lambda4[j] - (h/6)*(k14 + 2*k24 + 2*k34 + k44)
  lambda5[j-1]  <-  lambda5[j] - (h/6)*(k15 + 2*k25 + 2*k35 + k45)
  lambda6[j-1]  <-  lambda6[j] - (h/6)*(k16 + 2*k26 + 2*k36 + k46)
}

# Improved treatment function based upon optimality condition
temp <- (1/2)*((lambda2-lambda1)*EB/(SB+EB+IB) +(lambda3-lambda1)*IB/(SB+EB+IB))
u1 <- pmin(umax,pmax(0,temp,na.rm=TRUE),na.rm=TRUE)
u <- u1*(1-(.99)^(iter))+oldu*(0.99)^(iter)

# Checking for convergence
temp1 <- delta*sum(abs(u)) - sum(abs(oldu - u))
temp2 <- delta*sum(abs(SB)) - sum(abs(oldSB - SB))
temp3 <- delta*sum(abs(EB)) - sum(abs(oldEB - EB))
temp4 <- delta*sum(abs(IB)) - sum(abs(oldIB - IB))
temp5 <- delta*sum(abs(SS)) - sum(abs(oldSS - SS))
temp6 <- delta*sum(abs(ES)) - sum(abs(oldES - ES))
temp7 <- delta*sum(abs(IS)) - sum(abs(oldIS - IS))
temp8 <- delta*sum(abs(lambda1)) - sum(abs(oldlambda1 - lambda1))
temp9 <- delta*sum(abs(lambda2)) - sum(abs(oldlambda2 - lambda2))
temp10 <- delta*sum(abs(lambda3)) - sum(abs(oldlambda3 - lambda3))
temp11 <- delta*sum(abs(lambda4)) - sum(abs(oldlambda4 - lambda4))
temp12 <- delta*sum(abs(lambda5)) - sum(abs(oldlambda5 - lambda5))
temp13 <- delta*sum(abs(lambda6)) - sum(abs(oldlambda6 - lambda6))

test <- min(temp1, min(temp2, min(temp3, min(temp4, min(temp5,
        min(temp6, min(temp7, min(temp8,min(temp9, min(temp10, 
          temp11,min(temp12, temp13),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE),
          na.rm=TRUE),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE),na.rm=TRUE)
}

#Done....asssigning output of optimal vectors

  y <- matrix(0,8,M+1)
  y[1,] <- t
  y[2,] <- SB
  y[3,] <- EB
  y[4,] <- IB
  y[5,] <- SS
  y[6,] <- ES
  y[7,] <- IS
  y[8,] <- u

  y.df <- as.data.frame(t(y))
  colnames(y.df) <- c("time", "SB", "EB", "IB", "SS", "ES", "IS", "u")
  
  return(y.df)
  
}) # ends with(as.list...

} 

