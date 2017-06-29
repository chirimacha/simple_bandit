
'#####################################################
Running simulations for Global vs. Bandit


Creates and subdivides infestations and compares running the bandit vs. doing global ring search
Bandit search: bandit selects an arm and ring search (battleship) conducted in that arm with results reported 
to the bandit (and arm pulled) every X times 
Global search: Conduct ring search on a infestation created by combining the arms of the bandit into one large infestation

Calls "FunctionsForSimulatedSearch" which contains:
1) call to the infestation generator (create simulated data)
2) the code to conduct a ring search on simulated data
3) call to bandit parameters
4) Code to run simulated searches
'
#https://cran.r-project.org/web/packages/neldermead/neldermead.pdf


#TO DO: create an environmental variable that points to your simple_bandit repository
setwd(Sys.getenv("SIMPLE_BANDIT")) 

#Call simulated search function
source("code/FunctionsForSimulatedSearch.R")

#Infestation and search paramaters [per arm of the bandit]
'
@param max_cost:          limit for number of searches by RingSearch based on a cost. 
@param ncol:              ncol in infestation
@param nrow:              nrow in infestation
@param num_days:          number of days of infestation if not using prevalence 
@param num_starts:        number of foci that begin each infestation [FIXME: Set to number of potential reinfested foci]
@param p_hop:             probability of hop
@param p_skip:            probability of skip
@param p_jump:            probability of jump
@param pMovePerDay:       probability of movement per site
@param pClearancePerday: spontaneous clearance per day
@param params.bandit      c(learning,discount)
}
'

params_grid_sim <- list(  #based on latest from Barbu (4/14/2016)
  max_cost=10000,
  ncol=100,
  nrow=100,
  num_starts=30,
  p_hop=0.68,
  p_skip=0.16,
  p_jump=0.15,
  pMovePerDay=0.01*7.0, #per site.  to accel the sim, we use faster move rate than in the Barbu et al., but shorter time to complete the colonization (365 days)
  pClearancePerDay=0.000,  #spontaneous disappearance, per day
  ring=1,
  verbose=F
)

#Set parameters for prevalence in arms fo the bandit
#TODO: could replace with a simple vector, instead of named values
params.arm.TwoTwoTwo <- list (
  prev.arm1=.02,
  prev.arm2=.02,
  prev.arm3=.02
)

params.arm.TwoTwoTen <- list (
  prev.arm1=.02,
  prev.arm2=.02,
  prev.arm3=.10
)

params.arm.ZeroTen <- list (
  prev.arm1=.00,
  prev.arm2=.10
)


#Set the parameters of the Bandit

#BANDIT TYPE RC
##' RC Bandit internals 
#   @param learning: bandit learning rate
#   @param discount: bandit discount rate
# 
#   Bandit is updated with
#   F=preferences
#   r=reward
#   R=mean reward
#   A=discount factor
#   B=learning rate
# 
#   F(t+1)=F(t)+B(r(t)-R(t)) <- update first
#   R(t+1)=A(r(t))+(1-A)(R(t))
# 
params.bandit.rc.default <- list(
  name = "rc",
  discount_factor = 0.05,
  learning_rate   = 0.3444444
)
#TODO: are these the correct parameters?

#UCB bandit
params.bandit.ucb <- list(
  name = "ucb1"
)

max.cores.needed = 3  #needs to be set for each machine

##############################################################
print("STARTING comparisons ...")
##############################################################
# 
#Example Experiment
#ZeroTen vs. Bandit - power of the bandit with many starts that favor the bandit

results1 <- ZCompareGlobalVsBandit(params.bandit=params.bandit.rc.default, bandit.time=400, bandit.block.size=10,
                                   params.grid = params_grid_sim, params.arm=params.arm.ZeroTen,
                                   num.replications = 2, SearchStrategy = RingSearch,  altAlgo = RingSearch,
                                   cpu.cores.needed = max.cores.needed)


p1 <-  ZGraphResults.Formated(data=results1, params.grid=params_grid_sim, params.arm=params.arm.ZeroTen)
