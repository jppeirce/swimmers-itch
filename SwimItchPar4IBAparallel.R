## Rcode to accompany the paper
## A Mathematical Model for the Control of Swimmer's Itch
## J.P. Peirce^{1,3}, J.J. Pellett^1, G.J.  Sandland^2,3
## Code written by J. Peirce (jpeirce@uwlax.edu)
#####

library(deSolve)
library(foreach)
library(doParallel)
library(lhs)

### I. Mathematical Model
# Model Parameters
m.par <- c(
  rho = 0.745,  # prevalence of parasite in migrating merganser
  s = 0.40, # survival probability of egg to adult merganser
  beta = 2.13696*10^(-8), # transmission rate from infected snail to susceptible merganser
  muB = 1/(10*364.25), # natural death rate of merganser
  kB = 1.908*10^(-6), # mortality rate in merganser due to infection
  chi = 0.000035046144, # transmission rate from infected merganser to susceptible snail
  muS = (1-.333)/77, # natural death rate of snail
  gammaB = 1/14, # rate out of exposed bird population
  gammaS = 1/35, # rate out of exposed snail population
  kS = 5.48*10^(-3), # snail death rate due to parasite
  Tf = 7*30  # length of season
  )

# Initial Conditions
m.ic <- c(
  SB = 0,
  EB = 0,
  IB = 0,
  SS = 1000000,
  ES = 0,
  IS = 0
)

# Bird Migration Function
GammaIn <- function(t){
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

SwimItchModel <-function(t, state, parameters){
  with(as.list(c(state,parameters)), {
    dSB <- (1-rho)*GammaIn(t) + s*bB(t) - beta*SB*IS - muB*SB
    dEB <- beta*SB*IS - (gammaB+muB)*EB
    dIB <- rho*GammaIn(t) + gammaB*EB - (muB+kB)*IB
    dSS <- bS(t)*SS - chi*SS*IB - muS*SS
    dES <- chi*SS*IB-(gammaS+muS)*ES
    dIS <- gammaS*ES-(muS+kS)*IS
  list(c(dSB,dEB,dIB,dSS,dES,dIS))
  }) # ends with(as.list())
}

# To graph solutions:
dt <- 0.1 # timestep for ode solver
soln <- ode(y=m.ic, times=seq(0,m.par["Tf"], by=dt), func=SwimItchModel, parms=m.par)

### II. Parameter Estimation for chi, beta, and kB

# Read in data
snail.data <- read.csv('./PrevSnailData.csv',header=TRUE,sep=",")
names(snail.data) <- c("location", "date", "time", "prev.snail.mu", "prev.snail.SE")

higgins.snail.data <- subset(snail.data, subset = snail.data$location == "Higgins")
time.data <- higgins.snail.data$time
prev.snail.data <- higgins.snail.data$prev.snail.mu

# "Constant" Model Parameters
c.par <- c(
  rho = 0.745,  # prevalence of parasite in migrating merganser
  s = 0.40, # survival probability of egg to adult merganser
  muB = 1/(10*364.25), # natural death rate of merganser
  muS = (1-.333)/77, # natural death rate of snail
  gammaB = 1/14, # rate out of exposed bird population
  gammaS = 1/35, # rate out of exposed snail population
  kS = 5.48*10^(-3), # snail death rate due to parasite
  Tf = 7*30  # length of season
)
# Parameters to be estimated from data
# Initial guess
opt.par <-c(
  beta = 1*10^(-7), # transmission rate from infected snail to susceptible merganser
  chi = 1*10^(-5),
  kB = 1*10^(-6) # mortality rate in merganser due to infection
)
m.par <- c(c.par,opt.par) # both parameter sets together

# Least Squares function - sum of the square distances from model output to data
LeastSq <- function(opt.par){
  ssr <- NA
  if(all(opt.par > 0)){ # requires all parameters to remain not too small
    soln <- ode(y = m.ic, 
                times = seq(0,m.par["Tf"], by=dt),
                func=SwimItchModel, 
                parms=c(c.par,opt.par))
    t.de <- time.data*(1/dt)+1 # rescale t.data into de time frame
    prev.snail.sim <- soln[t.de,"IS"]/(soln[t.de,"SS"]+soln[t.de,"ES"]+soln[t.de,"IS"])
    ssr <- sum((prev.snail.data-prev.snail.sim)^2)
  }
  return(ssr)
}

### III. Genetic Algorithm Approach
#Multiple runs
for (i in 1:30){

start_time <- Sys.time()

# GA parameters
gen.total <- 5000 # number of individuals in each generation
percent.crossover <- 0.25
cross.total <- round(percent.crossover*gen.total)
percent.migration <- 0.10
mig.total <- round(percent.migration*gen.total)
mig.frequency <- 5 # number of years between migration events

### GA functions ###
## To mutate a single individual
mutator <- function(individual){
  # Identify the chromosomes (i.e. parameters) to mutate 
  mutate <- sample(x = 0:1, size = 3, replace = TRUE)
  # Comute a new set of chromosomes
  mut <- data.frame(beta = sample(par.matrix$beta, 1 ), 
                    chi = sample(par.matrix$chi, 1),
                    kB = sample(par.matrix$kB, 1))
  # Mutate (replace) those chomosomes which mutate is TRUE
  individual[mutate == TRUE] <- mut[mutate==TRUE]
  individual[4] <- LeastSq(individual[1:3])
  return(individual)
}

## Crossover Function - crosses two parameter trios 
offspring.total <- 2 # offspring produced per parent
crossover <- function(parent1, parent2, mutate = FALSE){
  # Create blank data frame
  offspring <- data.frame(beta = integer(offspring.total),
                          chi = integer(offspring.total),
                          kB = integer(offspring.total),
                          LS = integer(offspring.total))
  for (i in 1:offspring.total){
    # decide which parent contributes gene
    k <- sample(x = 1:2, size = 3, replace = TRUE)
    offspring[i,1:3] <- as.data.frame(ifelse(k==1, yes = parent1[1:3], no = parent2[1:3]))
    # Mutate and update this offspring's details if mutate = TRUE
    if (mutate) offspring[i,] <- mutator(individual = offspring[i,])
    offspring$LS[i] <- LeastSq(offspring[i,1:3])
  }
  return( offspring )
}
#############

## Initial Generation selected by Latin hypercube design
# Initial guess for parameters
N0 <- 100000
opt.par.init <- data.frame(
  beta = 1*10^(-7), # transmission rate from infected snail to susceptible merganser
  chi = 1*10^(-5), # transmission rate from infected merganser to susceptible snail
  kB = 1*10^(-6) # snail death rate due to parasite
)

# latin hyper cube random values 0 to 1
par.matrix <- as.data.frame(randomLHS(N0,3)) 
colnames(par.matrix) <- c("beta", "chi", "kB")

# rescale so that value are between 10^(-1) and 10^1
# times parameter values to create a large matrix of possible values
par.matrix$beta <- qunif(par.matrix$beta, 
                         min =  10^(-1)*opt.par.init$beta,
                         max =  10^(1)*opt.par.init$beta)
par.matrix$chi <- qunif(par.matrix$chi, 
                        min =  10^(-1)*opt.par.init$chi,
                        max =  10^(1)*opt.par.init$chi)
par.matrix$kB <- qunif(par.matrix$kB, 
                       min =  10^(-1)*opt.par.init$kB,
                       max =  10^(1)*opt.par.init$kB)

#select (randomly) the first generation)
gen <- par.matrix[sample(1:N0,gen.total, replace = FALSE),]
rownames(gen) <- 1:gen.total

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload computer
registerDoParallel(cl)

# compute the fitness (LS) of the first generation
gen$LS <- foreach(i=1:gen.total, .combine='c') %dopar% {
  library(deSolve)
  LeastSq(gen[i,1:3])
}

#stopCluster(cl) # closes the parallel cluster of cores

newOrder <- with(data = gen, order(LS)) # sorts the individuals by LS
gen.now <- gen[newOrder,]; rownames(gen.now) <- 1:gen.total

###### Now repeat for multiple generations
gen.max <- 100 # max number of generations
tol <- 10^(-8) # repeat while generational difference in mean(LS) > tol

diff.meanLS <- 1
generation <- 1
while ((generation < gen.max) & (diff.meanLS > tol)){
  generation <- generation + 1
  meanLS.now <- mean(gen.now$LS)
  # migration event
  if (generation%%mig.frequency ==0){
    mig <- data.frame( beta = sample(par.matrix$beta, mig.total, replace = FALSE), 
                       chi = sample(par.matrix$chi, mig.total, replace = FALSE),
                       kB = sample(par.matrix$kB, mig.total, replace = FALSE),
                       LS = rep(NA, length.out = mig.total))
    mig$LS <- foreach(i=1:mig.total, .combine='c') %dopar% {
      library(deSolve)
      LeastSq(mig[i,1:3])
    }
    gen.now <- rbind (gen.now, mig)
  }
  temp <- foreach(i=1:cross.total, .combine= rbind) %dopar% {
    library(deSolve)
    crossover(parent1 = gen.now[i,],
              parent2 = gen.now[sample( (1:gen.total)[-i], 1),], #random gene that is not ith
              mutate = TRUE)
  }
  gen.next <- rbind(gen.now,temp)
  newOrder <- order(gen.next[,"LS"])
  gen.next <- gen.next[newOrder,]; rownames(gen.next) <- 1:nrow(gen.next)
  gen.now <- gen.next[1:gen.total,]
  meanLS.next <- mean(gen.now$LS)
  diff.meanLS <- abs(meanLS.now-meanLS.next) 
}

stopCluster(cl) # closes the parallel cluster of cores

end_time <- Sys.time()
time <- end_time - start_time
run.output.df <- cbind.data.frame(head(gen.now,1),time,generation)

if (i==1){output.df <- run.output.df} 
  else {output.df <- rbind.data.frame(output.df,run.output.df)}

} #end multiple runs

write.csv(output.df, file = "optparam.csv")