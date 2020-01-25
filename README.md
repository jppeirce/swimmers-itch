# swimmers-itch
R code to accompany the paper, "A Mathematical Model for the Control of Swimmer's Itch."

File summary:
swimitchcontrol.R - Start here! This is the main file.  It contains the commands for all graphs included in our paper.
swimitchfun.R - Yes, this file is fun!! It contains the Forward-backward sweep method with a 4th order Runge-Kutta scheme to solve the    		state equations and their corresponding adjoint equations. 
SwimItchPar4IBAparallel.R - Code that applies a genetic algorithm scheme to the minimization process used to approximate the unknown 					parameters $\beta$, $\chi$, and $k_{\rm B}$. Designed for parallel computation.
optparam.csv - The output of SwimItchPar4IBAparallel.R contains the results of 30 trials. Contains best (lowest fitness score) parameter 			values, the computational time (in hours), and the number of generations before reaching tolerance.
PractOptCon4IBAparallel.R - Code that determines the day(s) treatment should be given when $M$ birds are treated.  Designed for parallel 			computation.
opttime1.csv and opttime2.csv - The output of PractOptCon4IBAparallel.R.
