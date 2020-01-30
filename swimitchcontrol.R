## Rcode to accompany the paper
## A Mathematical Model for the Control of Swimmer's Itch
## J.P. Peirce^{1,3}, J.J. Pellett^1, G.J.  Sandland^2,3
## Code written by J. Peirce (jpeirce@uwlax.edu)
#####

# NOTE: Remember to set working directory to current location
source("swimitchfun.R")
library(tidyverse)
library(latex2exp)
library(deSolve)
library(matrixStats)

# Optimal parameters found from genetic algorith approach (SwimItchPar4IBAparallel.R)
opt.par.sim <- read.csv('./optparam.csv',header=TRUE,sep=",")
m.par <- c(
  rho = 0.745,  # prevalence of parasite in migrating merganser
  s = 0.40, # survival probability of egg to adult merganser
  beta = mean(opt.par.sim$beta), # transmission rate from infected snail to susceptible merganser
  muB = 1/(10*364.25), # natural death rate of merganser
  kB = mean(opt.par.sim$kB), # mortality rate in merganser due to infection
  chi = mean(opt.par.sim$chi), # transmission rate from infected merganser to susceptible snail
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

## Figure 3 Functions that increase the bird population attributed to migration rate in(t)
## (birds/day) and hatchling survival rate sbB(t) (birds/day).
figure3.df <- data.frame(time = 1:m.par["Tf"], 
                         LambdaIn = LambdaIn(1:m.par["Tf"]),
                         sbB = m.par["s"]*bB(1:m.par["Tf"]))
ggplot(figure3.df) + 
  geom_line(aes(x=figure3.df$time, y=figure3.df$LambdaIn, linetype="\u039Bin")) + 
  geom_line(aes(x=figure3.df$time, y=figure3.df$sbB, linetype="s bB")) +
  # Legend order
  scale_linetype_manual(values = c("\u039Bin"='dashed',"s bB"='solid')) +
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  # Axis Formatting
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,6), breaks = seq(0,6,1), expand = c(0,0)) +
  # Labels 
  labs(x="Time (days)",
       y="Birds per day",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.8, .8),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=16),
        aspect.ratio = .7)
# Initial Conditions
m.ic <- c(
  SB = 0,
  EB = 0,
  IB = 0,
  SS = 1000000,
  ES = 0,
  IS = 0
)

### Solve the system of equations WITHOUT treatment
y.df  <-  SwimItchControl(m.par, m.ic, A = 0, umax = 0)

# Figure 4 A simulation of the parasite prevalence in the bird population during
#the summer residency season in untreated years

y.df$prev.bird <- y.df$IB/(y.df$SB+y.df$EB+y.df$IB)
y.df$prev.bird[1] <- m.par["rho"]
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure4.pdf", width = ((1+sqrt(5))/2)*7)
# Figure 4 A simulation of the parasite prevalence in the bird population during the summer 
#residency season in untreated years. Labeled times $tm$, $\tau$, and $\tau+tm$ correspond to 
#end of migration and the start and end of egg hatching period.
ggplot(data = y.df) + 
  geom_line(aes(x=time, y=prev.bird, linetype="Sim. Prev. in Bird"),show.legend = FALSE) + 
  geom_abline(slope=0, intercept = m.par["rho"], linetype="dashed",show.legend = FALSE)+
  geom_abline(slope=0, intercept = tail(y.df$IB/(y.df$SB+y.df$EB+y.df$IB),1), linetype="dashed",show.legend = FALSE)+
  geom_segment(aes(x=60, y=0, xend=165, yend=0), size=1.5, show.legend = FALSE) +
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.1), expand = c(0,0)) +
  # Labels 
  labs(x="Time (days)",
       y="Sim. Prev. in Birds",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.8, .9),
        legend.key.width = unit(1.5,"cm"),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=16),
        aspect.ratio = .7)
#dev.off()

# Read in Higgin's Lake data
snail.data <- read.csv('./PrevSnailData.csv',header=TRUE,sep=",")
names(snail.data) <- c("location", "date", "time", "prev.snail.mu", "prev.snail.SE")
higgins.snail.data <- subset(snail.data, subset = snail.data$location == "Higgins")
time.data <- higgins.snail.data$time
prev.snail.data <- higgins.snail.data$prev.snail.mu

# Figure 5: (a) Simulated parasite prevalence in the snail population without 
# treatment using infection data from Higgins Lake, MI. Data includes mean (and one 
# standard error) prevalence across 10 sites at Higgins Lake at four summer time points. 
# Horizontal dashed lines are at ideal and epidemic snail prevalence thresholds. 
# (b) Simulated merganser densities in a season without treatment.
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure5a.pdf",width = ((1+sqrt(5))/2)*7)
ggplot(data = y.df) + 
  # Plot population density graphs 
  geom_line(aes(x=time, y=IS/(SS+ES+IS), linetype="Sim. Prev. in Snail"),show.legend = FALSE) + 
  geom_point(data = higgins.snail.data, 
             aes(x = time.data, y = prev.snail.data)) + 
  geom_errorbar(data = higgins.snail.data,
                aes(x = time.data,
                    ymin = higgins.snail.data$prev.snail.mu - higgins.snail.data$prev.snail.SE, 
                    ymax = higgins.snail.data$prev.snail.mu + higgins.snail.data$prev.snail.SE)) +
  geom_segment(aes(x=60, y=0, xend=165, yend=0), size=1.5, show.legend = FALSE) +
  geom_abline(slope=0, intercept = 0.02, linetype="1F",show.legend = FALSE)+
  geom_abline(slope=0, intercept = 0.0024, linetype="1F",show.legend = FALSE)+
   # Legend order
  scale_shape_manual(values = c("Sim. Prev. in Snail"='solid')) +
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  # Axis Formatting
  scale_x_continuous(limits = c(0,210), breaks = seq(0,m.par["Tf"],50), expand = c(0,0))+
  scale_y_continuous(limits = c(0,.1), breaks = seq(0,.1,.010), expand = c(0,0)) +
  # Labels 
  labs(x="Time (in days)",
       y="Prevalence in Snail Host",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.3, .9),
        legend.key.width = unit(1.5,"cm"),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure5b.pdf",width = ((1+sqrt(5))/2)*7)
ggplot(data = y.df,
       aes( x = time, y = SB, linetype = "SB")) + 
  geom_line() + 
  geom_line( aes( x = time, y = IB, linetype = "IB")) +
  geom_line( aes( x = time, y = EB, linetype = "EB")) +
  scale_linetype_manual(breaks = c("SB","EB","IB"),
                        values = c(SB = "solid", EB = "dotted", IB = "dashed")) + 
  geom_segment(aes(x=60, y=0, xend=165, yend=0), size=1.5, show.legend = FALSE) +
  labs(color = "Legend",
       x = "Time (in days)",
       y = "Bird Density",
       linetype = "Legend") +
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,90), breaks = seq(0,90,20), expand = c(0,0))+
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  theme_classic()+
  theme(legend.position = c(.9, .5),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()

### Sensitivity during resident period
# Residency period
res.start <- 60
res.end <- 165
meanIS0 <- mean(y.df$IS[which(y.df$time >= res.start & y.df$time <= res.end)])
meanISnew <- data.frame(
  up = numeric(length = (length(m.par)-1)),
  down = numeric(length = (length(m.par)-1))
)
m.par.up <- m.par
m.par.down <- m.par
# Computes the mean IS during residency period when each parameter is increase by 10%
for (i in 1:(length(m.par)-1)){ # not computing index for Tf
  m.par.up[i] <-  1.10 * m.par[i]
  m.par.down[i] <-  0.90 * m.par[i]
  y.df.up  <-  SwimItchControl(m.par.up, m.ic, A = 0, umax = 0)
  y.df.down  <-  SwimItchControl(m.par.down, m.ic, A = 0, umax = 0)
  meanISnew$up[i] <- mean(y.df.up$IS[which(y.df$time >= res.start & y.df$time <= res.end)])
  meanISnew$down[i] <- mean(y.df.down$IS[which(y.df$time >= res.start & y.df$time <= res.end)])
  m.par.up <- m.par
  m.par.down <- m.par
}
# without latex par = (names(m.par)[1:(length(m.par)-1)])
# Use unicode for greek letters
sensitivity.index <- data.frame(
  par =c("\u03C1","s","\u03B2","\u03BCB","kB","\u03A7","\u03BCS",
         "\u03B3B","\u03B3S","kS"),
  index = 100*colMaxs(rbind(abs((meanISnew$up-meanIS0)/meanIS0),abs((meanISnew$down-meanIS0)/meanIS0))),
  Host = c("bird","bird","parasite","bird",
           "bird","parasite","snail","bird","snail","snail")
)
# Sort   
sensitivity.index <- sensitivity.index[order(-sensitivity.index$index),]
ggplot(data = sensitivity.index, 
       aes(x=reorder(par,-index), y=index, fill = Host)) +
  geom_bar(stat = "identity") +
  scale_fill_grey() +
  scale_y_continuous(limits = c(0,14), breaks = seq(0,14,1), expand = c(0,0))+
  labs(x="Model Parameters",
       y="Sensitivity Index",
       linetype = "Legend") +
  theme_classic() +
  theme(legend.position = c(.8, .8),
        #    legend.key.width = unit(1.5,"cm"),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=16),
        aspect.ratio = .7)

rm(y.df.up, y.df.down)

### Optimal Control
# Solve optimal control problem with the two weight factors A
umax <- round(max(y.df$SB+y.df$EB+y.df$IB)) # max birds
Alo <- 1.23*10^(-5)  # low weight factor
ylo.df  <-  SwimItchControl(m.par, m.ic, Alo, umax)
Ahi <- 45.75*10^(-5) # high weight factor
yhi.df  <-  SwimItchControl(m.par, m.ic, Ahi, umax)

# Check prevalence in snail is below thresholds
100*ylo.df$IS[floor(5.5/7*length(ylo.df$time))]/(ylo.df$SS[floor(5.5/7*length(ylo.df$time))]
                                                 +ylo.df$ES[floor(5.5/7*length(ylo.df$time))]
                                                 +ylo.df$IS[floor(5.5/7*length(ylo.df$time))])
100*yhi.df$IS[floor(5.5/7*length(yhi.df$time))]/(yhi.df$SS[floor(5.5/7*length(ylo.df$time))]
                                                 +yhi.df$ES[floor(5.5/7*length(ylo.df$time))]
                                                 +yhi.df$IS[floor(5.5/7*length(ylo.df$time))])

# Figure 7 (a) Optimal control $u(t)$ with weights $A$. (b) Simulated density 
# of infected bird in the three cases of treatment. (c) Simulated prevalence of 
# the snail density in the three cases of treatment. Horizontal dashed lines at
#ideal and epidemic snail prevalence thresholds.
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure7a.pdf",width = ((1+sqrt(5))/2)*7)
ggplot(data = ylo.df) + 
  geom_line(aes(x=time, y = u, linetype="Alo")) +
  geom_line(data = yhi.df, aes(x=time, y = u, linetype="Ahi")) +
  scale_linetype_manual(values = c(Alo='dotted',Ahi='dashed')) +
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  # Axis Formatting
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
    scale_y_continuous(limits = c(0,4), breaks = seq(0,4,.5), expand = c(0,0)) +
  # Labels 
  labs(x="Time (in days)",
       y="Treatment",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.8, .8),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()

#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure7b.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = y.df) + 
  geom_line(aes(x=time, y = IB, linetype="No Treatment")) +
  geom_line(data = yhi.df, aes(x=time, y = IB, linetype="Ahi")) +
  geom_line(data = ylo.df, aes(x=time, y = IB, linetype="Alo")) +
  scale_linetype_manual(values = c("No Treatment"='solid', Alo='dotted',Ahi='dashed')) +
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  # Axis Formatting
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
    scale_y_continuous(limits = c(0,45), breaks = seq(0,45,10), expand = c(0,0)) +
  # Labels 
  labs(x="Time (in days)",
       y="Infected Bird",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.3, .8),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure7c.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = y.df) + 
  geom_line(aes(x=time, y = IS/(SS+ES+IS), linetype="No Treatment")) +
  geom_line(data = ylo.df, aes(x=time, y = IS/(SS+ES+IS), linetype="Alo")) +
  geom_line(data = yhi.df, aes(x=time, y = IS/(SS+ES+IS), linetype="Ahi")) +
  geom_point(data = higgins.snail.data, 
             aes(x = time.data, y = prev.snail.data)) + 
  geom_errorbar(data = higgins.snail.data,
             aes(x = time.data,
                 ymin = higgins.snail.data$prev.snail.mu - higgins.snail.data$prev.snail.SE, 
                ymax = higgins.snail.data$prev.snail.mu + higgins.snail.data$prev.snail.SE)) +
  geom_segment(aes(x=60, y=0, xend=165, yend=0), size=2, show.legend = FALSE) +
  scale_linetype_manual(values = c("No Treatment"='solid', Alo='dotted',Ahi='dashed')) +
  geom_abline(slope=0, intercept = 0.02, linetype="1F",show.legend = FALSE)+
  geom_abline(slope=0, intercept = 0.0024, linetype="1F",show.legend = FALSE)+
  guides(linetype=guide_legend(keywidth = 3, keyheight = 1), color=FALSE) +
  # Axis Formatting
  scale_x_continuous(limits = c(0,210), breaks = seq(0,m.par["Tf"],50), expand = c(0,0))+
  scale_y_continuous(limits = c(0,.1), breaks = seq(0,.1,.010), expand = c(0,0)) +
  # Labels 
  labs(x="Time (in days)",
       y="Prevalence in Snail Host",
       linetype = "Legend") + 
  # Plotting Theme
  theme_classic() +
  theme(legend.position = c(.3, .8),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()

###### Longevity of Optimal Treatment
### To graph pi (prevalence after one year) vs rho0 (initial prevalance)
rho0N <- 10 # number of rho0 to initially use between observed prevalence range
rho0vec <- seq( from = 0.6, to = 0.89, length.out = rho0N )

pivec <- matrix(0,1,rho0N)
for (i in 1:rho0N){
  m.par["rho"] <- rho0vec[i]
  y.df  <-  SwimItchControl(m.par, m.ic, A = 0, umax = 0)
  pivec[i] <- tail(y.df$EB+y.df$IB,n=1)/(tail(y.df$SB,n=1)+tail(y.df$EB,n=1)+tail(y.df$IB,n=1));
}
# plot(rho0vec,pivec,type="p", main="rho0 vs. pi")

# Create a dataframe whose first column is rho0 and second column is pi
pifunc <- data.frame(rho0=rho0vec, pi=t(pivec))
# Create a linear regression model
pifunc.lm <- lm(pi~rho0,data=pifunc)
# Display coefficients of model
pifunc.lm$coefficients
summary(pifunc.lm)
# Plot simulated data with linear model
#plot(pifunc,type="p", main="rho0 vs. pi")
#abline(pifunc.lm$coefficients)

# prev in bird population at end of season
tail((y.df$EB+y.df$IB)/(y.df$SB+y.df$EB+y.df$IB),1)
tail((ylo.df$EB+ylo.df$IB)/(ylo.df$SB+ylo.df$EB+ylo.df$IB),1)
tail((yhi.df$EB+yhi.df$IB)/(yhi.df$SB+yhi.df$EB+yhi.df$IB),1)

# Longevity of Treatment - ALPHA 0.80
d0 <- pifunc.lm$coefficients[2] # slope of linear model when there is no treatment
alpha <- 0.80 # percent of original prevalence
rho <- 0.745
pi0 <- d0*rho + pifunc.lm$coefficients[1]

piAlo <- tail(ylo.df$EB+ylo.df$IB,n=1)/(tail(ylo.df$SB,n=1)+tail(ylo.df$EB,n=1)+tail(ylo.df$IB,n=1))
piAhi <- tail(yhi.df$EB+yhi.df$IB,n=1)/(tail(yhi.df$SB,n=1)+tail(yhi.df$EB,n=1)+tail(yhi.df$IB,n=1))
(TreatTimelo <- 1+(1/log(d0))*log(((1-alpha)*rho)/(pi0-piAlo)))
(TreatTimehi <- 1+(1/log(d0))*log(((1-alpha)*rho)/(pi0-piAhi)))

# Longevity of Treatment - ALPHA = 0.95
d0 <- pifunc.lm$coefficients[2] # slope of linear model when there is no treatment
alpha <- 0.95 # percent of original prevalence
rho <- 0.745
pi0 <- d0*rho + pifunc.lm$coefficients[1]

piAlo <- tail(ylo.df$EB+ylo.df$IB,n=1)/(tail(ylo.df$SB,n=1)+tail(ylo.df$EB,n=1)+tail(ylo.df$IB,n=1))
piAhi <- tail(ylo.df$EB+yhi.df$IB,n=1)/(tail(yhi.df$SB,n=1)+tail(yhi.df$EB,n=1)+tail(yhi.df$IB,n=1))
(TreatTimelo <- 1+(1/log(d0))*log(((1-alpha)*rho)/(pi0-piAlo)))
(TreatTimehi <- 1+(1/log(d0))*log(((1-alpha)*rho)/(pi0-piAhi)))

###### Practical Optimal Control
# Differential equation model of equations without control function

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
  if (treat1 == 0){
    return(as.data.frame(ode(y=m.ic, times=seq(0,m.par["Tf"], by=dt), func=SwimItchModel, parms=m.par)))
  } else {
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
  }
  m.ic <- unlist(m.ic, row)
  #solve after treatment
  soln.afterT <- as.data.frame(
    ode(y=m.ic, times=seq(treat1, m.par["Tf"], by=dt), 
        func=SwimItchModel, parms=m.par)
  )
  return(rbind(soln.beforeT[soln.beforeT$time<treat1,], soln.afterT))
}

PracticalModel2 <- function(treat1, treat2, max.treat){
  #solve until treatment
  soln.beforeT1 <- as.data.frame(
    ode(y=m.ic, times=seq(0,treat1, by=dt), func=SwimItchModel, parms=m.par)
  )
  # numbers of at treatment day
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

#From parallel computation found in PractOptCon4IBAparallel.R
opt.time1 <- read.csv('./opttime1.csv',header=TRUE,sep=",")
opt.time2 <- read.csv('./opttime2.csv',header=TRUE,sep=",")
opt.time1 <- opt.time1[-1]
opt.time2 <- opt.time2[-1]

M.opt <- 22; t.opt <- opt.time1$day[opt.time1$M == M.opt] # for example

# FIGURE 8 (a) The practical optimal control of minimizing the 
# mean IS during recreational period by treating a maximum number of birds 
# $M$ on day $t_1$,  (b), (c) example of simulated infected bird and 
# infected snail treated(dashed) vs untreated (solid) when $M=22$ birds are 
# treated on day $t_1 = 13$.
notreat <- PracticalModel1(treat1 = 0, max.treat = NA)
NB <- (notreat$SB+notreat$EB+notreat$IB)
time2reach <- numeric(40)
for (m in 1:40){
  # compute time bird pop first reaches size m
  time2reach[m] <- as.numeric(notreat[which(NB > m)[1],]["time"])
}
time2reach.df <- data.frame(
  "day" = time2reach,
  "M" = c(1:40)
)
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure8a.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = opt.time1[1:40,],
       aes( x = day, y = M, linetype = "M" ),show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  geom_abline(slope=0, intercept = M.opt, linetype="1F",show.legend = FALSE)+
  geom_line(data = time2reach.df,
            aes( x= day, y = M, linetype = "min. time to reach M"),show.legend = FALSE)+
  scale_linetype_manual(values = c("M"='solid', "min. time to reach M" = 'dashed')) +
  labs(linetype = "Legend",
       x = "Day Treated",
       y = "Maximum Birds Treated") +
  scale_x_continuous(limits = c(0,30), breaks = seq(0,30,5), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,5), expand = c(0,0)) +
  theme_classic() +
  theme(legend.position = c(.3, .9),
        legend.key.width = unit(0.5,"cm"),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()

# # Example graph
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure8b.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = PracticalModel1(treat1 = t.opt, max.treat = M.opt),
       aes( x = time, y = IB, linetype = "IB treated")) +
  geom_line() +
  geom_line(data = y.df,
            aes( x = time, y = IB, linetype = "IB")) +
  scale_linetype_manual(breaks = c("IB treated","IB"),
                     values = c("IB treated" = 'dashed', IB = 'solid')) +
  labs(linetype = "Legend",
       x = "Time (in days)",
       y = "Bird Density") +
  scale_x_continuous(limits = c(0,210),
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
  #scale_x_continuous(limits = c(0,210), breaks = seq(0,m.par["Tf"],50), expand = c(0,0))+
  scale_y_continuous(limits = c(0,40), breaks = seq(0,40,10), expand = c(0,0))+
  theme_classic()+
  theme(legend.position = c(.8, .5),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()

#### treat 2
M.opt <- 15 # for example
t1.opt <- opt.time2$day1[opt.time2$M == M.opt]
t2.opt <- opt.time2$day2[opt.time2$M == M.opt]

#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure9a.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = opt.time2,
       aes( x = day1, y = M, shape = "day1", color = "day1" )) +
  geom_point(show.legend = FALSE) +
  geom_point(aes( x = day2, y = M, shape = "day2", color = "day2"),show.legend = FALSE)+
  scale_shape_manual( values = c(2,6) )+
  scale_color_manual( values= c('black', 'black'), guide = FALSE)+
  labs(shape = "Legend",
       x = "Day Treated",
       y = "Maximum Birds Treated") +
  geom_abline(slope=0, intercept = M.opt, linetype="dashed",show.legend = FALSE)+
  scale_x_continuous(limits = c(0,28), breaks = seq(0,28,2), expand = c(0,0))+
  scale_y_continuous(limits = c(0,43), breaks = seq(0,42,5), expand = c(0,0))+
  theme_classic() +
  theme(
#  legend.position = c(.2, .8),
#        legend.key.width = unit(1.5,"cm"),
#        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()
#pdf("~/OneDrive - University of Wisconsin-La Crosse/SwimmersItch/paper/figures/Figure9b.pdf", width = ((1+sqrt(5))/2)*7)
ggplot(data = PracticalModel2(treat1 = t1.opt, treat2 = t2.opt, max.treat = M.opt),
       aes( x = time, y = IB, linetype = "IB treated")) + 
  geom_line() + 
  geom_line(data = y.df,
            aes( x = time, y = IB, linetype = "IB")) +
  scale_linetype_manual(breaks = c("IB treated", "IB"),
                     values = c("IB treated" = 'dashed', IB = 'solid')) + 
  labs(linetype = "Legend",
       x = "Time (in days)",
       y = "Infected Bird Density") +
  scale_x_continuous(limits = c(0,210), 
                     breaks = c(21,31.5,52.5,210),
                     labels = c(TeX('$t_m$'),TeX('$\\tau$'),
                                TeX('$\\tau+t_m$'),"T"),
                     expand = c(0,0))+
  scale_y_continuous(limits = c(0,40), breaks = seq(0,50,10), expand = c(0,0))+
  theme_classic()+
  theme(legend.position = c(.8, .5),
        legend.key.width = unit(1.5,"cm"),
        legend.title = element_blank(),
        plot.margin = margin(10, 10,10, 10),
        text = element_text(size=18),
        aspect.ratio = .7)
#dev.off()