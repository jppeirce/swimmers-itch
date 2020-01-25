## R code to accompany the paper, "A Mathematical Model for the Control of Swimmer's Itch."

Swimmer's itch is an emerging disease caused by flatworm parasites that often use water birds as definitive hosts. When parasite larvae penetrate human skin they initiate localized inflammation that leads to intense itching.  Concerns about this issue have been growing recently due to an apparent increase in the global occurrence of swimmer's itch and its subsequent impacts on recreational activities and associated revenues.  Past work has identified the common merganser as a key definitive host for these worms in the United States; a number of snail species serve as intermediate hosts. Although previous attempts at controlling swimmer's itch have targeted snails, a handful of efforts have concentrated on treating water birds with the anthelmintic drug, praziquantel.  We construct a mathematical model of swimmer's itch and its treatment within the infected merganser population. Our goal is to identify merganser treatment regimes that minimize the number of infected snails thereby reducing the risk of human infections.  Optimal control of bird hosts is defined analytically and we include numerical simulations assuming different resource-allocation strategies. Results from the study may help to identify treatment protocols that lower merganser infection rates and ultimately reduce the occurrence of swimmer's itch in freshwater systems throughout the Midwest.


# File summary:

swimitchcontrol.R - Start here! This is the main file.  It contains the commands for all graphs included in our paper.

swimitchfun.R - Yes, this file is fun!! It contains the Forward-backward sweep method with a 4th order Runge-Kutta scheme to solve the    		state equations and their corresponding adjoint equations. 

SwimItchPar4IBAparallel.R - Code that applies a genetic algorithm scheme to the minimization process used to approximate the unknown 					parameters $\beta$, $\chi$, and $k_{\rm B}$. Designed for parallel computation.
optparam.csv - The output of SwimItchPar4IBAparallel.R contains the results of 30 trials. Contains best (lowest fitness score) parameter 			values, the computational time (in hours), and the number of generations before reaching tolerance.
PractOptCon4IBAparallel.R - Code that determines the day(s) treatment should be given when $M$ birds are treated.  Designed for parallel 			computation.
opttime1.csv and opttime2.csv - The output of PractOptCon4IBAparallel.R.
