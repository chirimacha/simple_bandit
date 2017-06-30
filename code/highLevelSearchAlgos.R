'
################################
Bandit Simulations to use to compare global and bandit searching

This code:
  -Calls Bandit Code 
  -[Optionally] Simulates infestations
  -Simulates Searches: Provides mechanism for both a random search and a ring search 
  -Simulates a search using a bandit based on the parameters of the bandit

To call this code you need to:

Set prevalence/arm parameters

######################################
'

library(reshape2)
library(ggplot2)
library(pROC)
library(sm)

setwd(Sys.getenv("SIMPLE_BANDIT"))

#CALL BANDIT CODE - This uses the bandit code created by S. Gutfraind and S. Nutman 
source("code/bandit.R")
source("code/helperFunctions.R")

##' Ring Search ("Battleship")
##'   1. Randomly selects sites until makes a hit
##'   2. Explores the ring around the known hit.  If exhausts all known rings, reverts to random search
##'   st = state.  on first call st==NULL which causes this to initialize using RingSearchInitialize
##'   random_search=TRUE puts the code into purely random search, i.e. does not use rings (for benchmarking)
##' @param infestation = a rectangular grid of sites
##' @param st = the state of the algorithm, including number of sites uncovered etc
##' @param max_actions = number of searches the algorithm can conduct, before exiting
##' @param random_search = visit at random, instead of going for rings (neighbors) of known infestations
##' @params params = parameters that describe size of grid, parameters for the infestation, and size of the ring to be searched

RingSearch <- function(infestation, st=NULL, max_actions=Inf, params=NULL, random_search=FALSE) {
  if(is.null(st)) {
    st <- RingSearchInitialize(infestation=infestation, params=params)
  }
  
  initial_cost <- tail(st$running_stats[["total_cost"]], 1)
  ContinueInspection <- params[["ContinueInspection"]]
  if(is.null(ContinueInspection)) {
    ContinueInspection <- function(infestation, latest, st) {
      return(latest$total_cost < params$max_cost & dim(latest)[1] < st$total_squares);
    }
  }
  next_stat <- tail(st$running_stats, 1)
  while(next_stat$total_visited < st$total_squares & next_stat$total_cost < initial_cost + max_actions & 
        ContinueInspection(infestation, next_stat, st)) {
    #next_stat$step <- next_stat$step + 1
    next_site <- NULL
    #if(next_suspected_nb < dim(suspected_nbs)[1]){
    #  browser()
    #}      
    while (st$next_suspected_nb <= dim(st$suspected_nbs)[1] & (! random_search)) {
      next_nb <- st$suspected_nbs[st$next_suspected_nb,]
      st$next_suspected_nb <- st$next_suspected_nb + 1
      if(st$visited[next_nb$lat, next_nb$lon] == 0) {
        next_site <- next_nb
        st$visited_coordinates<-rbind(st$visited_coordinates,next_site)
        break
      } else {
      }
    }
    while (is.null(next_site) & st$next_random_idx <= st$total_squares) {
      next_site <- st$randomized_sites[st$next_random_idx,]
      if (st$visited[next_site$lat, next_site$lon] == 0) {
        st$next_random_idx <- st$next_random_idx + 1
        st$visited_coordinates<-rbind(st$visited_coordinates,next_site)
        break
      } else {
        next_site <-NULL
      } 
      st$next_random_idx <- st$next_random_idx + 1
    }
    if(is.null(next_site)) {
      break
    }
    next_stat$total_cost    <- next_stat$total_cost + 1
    next_stat$total_visited <- next_stat$total_visited + 1
    st$visited[next_site$lat, next_site$lon] <- 1.0  #visitted == 1 (visited), == 2 (found infestation)
    if(infestation[next_site$lat, next_site$lon] > 0) {
      st$visited[next_site$lat, next_site$lon] <- 2.0 #dim(st$running_stat)[1]/10.0
      st$known_infested     <- rbind(st$known_infested, next_site)
      next_stat$total_found <- next_stat$total_found + 1
      next_stat$total_bugs  <- next_stat$total_bugs + (infestation[next_site$lat, next_site$lon])
      neighbors             <- RingSearchGetUnvisitedNbs(next_site, st$visited, ring=params$ring)
      st$suspected_nbs      <- rbind(st$suspected_nbs, neighbors)
      #debugg addition
      if (dim(st$suspected_nbs)[1]>0) {
        rownames(st$suspected_nbs) <- seq(dim(st$suspected_nbs)[1])
      }  
      # print(st$suspected_nbs)
      #print(sprintf("1, lon=%d, lat=%d", next_site$lon, next_site$lat))
    } else {
      #print(0)
    }
    #st$running_stats <- rbind(st$running_stats, next_stat)
    st$running_stats <- rbind(st$running_stats, next_stat)
    next_site <- NULL
  }
  st$running_stats$infestation_lower_bound <- st$running_stats$total_found/st$total_squares
  st$running_stats$infestation_estimate    <- st$running_stats$total_found/st$running_stats$total_visited
  row.names(st$running_stats) <- seq(dim(st$running_stats)[1])
  st$running_stats$unfoundprevalence <- (st$true.prevalence-st$running_stats$total_found)/(st$total_squares-st$running_stats$total_visited)
  
  if(params[["verbose"]]) {
    print(tail(st$running_stats,1))
  }
  return(st)
}

#Find unvisited houses in the ring around the infected house
RingSearchGetUnvisitedNbs <- function(house, visited,ring) {
  max_lat <- dim(visited)[1]
  max_lon <- dim(visited)[2]
  nbs <- read.csv(text="lon,lat")
  for(x in seq(house$lon-ring, house$lon+ring)) {
    if (x < 1 | x > max_lon) {
      next
    }
    for(y in seq(house$lat-ring, house$lat+ring)) {
      if (y < 1 | y > max_lat) {
        next
      }
      nb = list(lon=x,lat=y)  #wishlist: prioritize by range
      if (all(nb == house)) {
        next
      }
      if (visited[nb$lat,nb$lon] > 0) {
        next
      }
      nbs <- rbind(nbs, nb)
    }
  }
  return(nbs)
}

#Initialize parameters for ring search 
RingSearchInitialize <- function(infestation, params=NULL) {
  st <- list()
  st$total_squares <- dim(infestation)[1] * dim(infestation)[2]
  st$true.prevalence <- sum(colSums(infestation !=0))
  st$visited <- matrix(0, nrow=dim(infestation)[1], ncol=dim(infestation)[2])
  st$randomized_sites <- melt(st$visited)
  st$randomized_sites <- melt(st$visited)
  st$randomized_sites <- st$randomized_sites[sample(st$total_squares),]
  names(st$randomized_sites)<-c("lat", "lon", "val")
  st$randomized_sites$val <- NULL
  st$running_stats <- data.frame(total_cost=c(0),total_visited=c(0),total_found=c(0),total_bugs=c(0))
  
  st$known_infested <- read.csv(text="UNICODE")
  st$suspected_nbs  <- read.csv(text="UNICODE")
  st$next_suspected_nb <- 1
  st$next_random_idx <- 1    
  st$visited_coordinates <-read.csv(text="lon,lat")
  return(st)
}

##' Random Search 
##' Uses the RingSearch Algorithm to select a random search
RandomSearch <- function(infestation, st=NULL, max_actions=Inf, params=NULL, random_search=FALSE) {
  return(RingSearch(infestation=infestation, st=st, max_actions=max_actions, params=params, random_search=TRUE))
}



#HELPER FUNCTIONS FOR THE BANDIT

#Function to find reward for bandit
#Current: calculates rewards for a vector of results [number of searches]
#Rewards are log10(1+number of bugs found)
##' @params new_st: state of the bandit off of which rewards are being calculated
##' @params block.size: number of houses being searched for single pull of bandit
BugReward <- function(new_st, block.size) {
  last <- tail(new_st$running_stats$total_bugs, (block.size+1)) #using total_bugs instead of total_found 
  last.reward <- rep(0,block.size)
  for (k in 2:(block.size+1)) {
    j=k-1
    last.reward[j] <- log10(1+last[k]-last[j])
  }
  return(last.reward)
}

#Option could use houses found as reward
##' @params new_st: state of the bandit off of which rewards are being calculated
##' @params block.size: number of houses being searched for single pull of bandit

HouseReward <- function(new_st, block.size) {
  last <- tail(new_st$running_stats$total_found, (block.size+1)) #using total_found 
  last.reward <- rep(0,block.size)
  for (k in 2:(block.size+1)) {
    j=k-1
    last.reward[j] <- last[k]-last[j]
  }
  return(last.reward)
}

#Find total bugs: Determines the total bugs and total infested houses found over the block search
##' @params new_st: state of the bandit off of which rewards are being calculated
##' @params block.size: number of houses being searched for single pull of bandit

BugsFound <-function(new_st, block.size) {  
  last.bug      <- tail(new_st$running_stats$total_bugs,(block.size+1))
  last.house    <- tail(new_st$running_stats$total_found,(block.size+1))
  last.bug      <-last.bug[block.size+1]-last.bug[1]
  last.house    <-last.house[block.size+1]-last.house[1]
  bugs.found    <- (cbind(last.bug,last.house))
  return(bugs.found)
}


#Function to find remaining prevalence (can be used for testing)
#   UnfoundPrevalence <- function(new_st) {
#     unfound <-tail(new_st$running_stats$unfoundprevalence,1)
#     return(unfound)
#   }


##' Simulate a Bandit search on grid with blocks
##' @params test.time:        number of times to run the bandit
##' @params params:           infestation grid parameters
##' @params params.arm:       arms and prevalence parameters
##' @params infestation:      generate new infestations (if NULL), or start with old ones [for benchmarking]
##' each infestation is a rectangular grid
##' @params block.size: number of searches conducted each time an arm is pulled
##' @params SearchStrategy   strategy used to search (Ring/Random)
##' @params RewardFunction: function used to calculate rewards for the bandit (takes in new_st, block.size)
##' NOTE: this algorithm is designed to run with a UCB1 bandit or an RC bandit. If additional bandit algorithms are added,
##' addtional code will be needed to initialize that bandit in this function

BanditSearchGridBlocks <- function(test.time=NULL, params=NULL, 
                                   block.size=NULL, params.bandit=NULL, 
                                   params.arm=NULL, infestations=NULL,
                                   SearchStrategy=NULL,RewardFunction=NULL) {
  
  #Optionally set up infestations  
  if(is.null(infestations)) {
    #GENERATE INFESTATIONS
    infestations<- lapply(params.arm, GenerateInfestation,params=params) #new infestations for each arm of the bandit
  }
  n_arms <- length(infestations)
  search_stats <- list()
  
  #Information about the infestations generated (e.g. houses with bugs, total number of bugs on grid)
  infestation.stats <-InfestationStats(infestations =infestations) #See Helper functions for further documentation
  
  #SET UP BANDIT ARMS
  #initialize function for all infestations 
  search_stats <- lapply(infestations, SearchStrategy, st=NULL, max_actions=0, params=params)
  
  #Internal function to pull arm of the bandit
  pull_arm <- function(chosen_arm, search_stats,block.size, SearchStrategy) {
    st     <- search_stats[[chosen_arm]]
    new_st <- SearchStrategy(data.frame(infestations[[chosen_arm]]),st=st,max_actions=block.size,params=params) 
    search_stats[[chosen_arm]] <- new_st
    return(search_stats)
  }
  
  print("Stochastic search with bandit");
  
  #Setting up running tables to keep track of searches
  times       <- seq(1,test.time)
  arms        <- rep(0,length(times))
  rewards     <- rep(0,length(times))
  mR          <- rep(0,length(times))
  # unfoundprev <- rep(0,length(times))
  blocking <- rep(0,length(times))
  # ps_holding <- matrix(nrow=length(times),ncol=length(params.arm))
  bugs.houses <-matrix(nrow=length(times),ncol=2) #uncommented to track
  static.houses.infested <-rep(infestation.stats[,"total.houses.infested"],length(times))
  static.bugs.available <-rep(infestation.stats[,"total.bugs.available"],length(times))
  static.rewards.available <-rep(infestation.stats[,"total.rewards.available"],length(times))
  
  #Initialize the bandit 
  
  if (params.bandit$name =="rc") {
    bandit <- initialize_rc(n_arms=n_arms, learning_rate=params.bandit$learning_rate, discount_factor=params.bandit$discount_factor)
  }
  else if (params.bandit$name =="ucb1") {
   bandit <- initialize_ucb1(n_arms=n_arms, reward_func=params.bandit$reward_func)
  } else if (params.bandit$name=="egreedy") {
      bandit <-initialize_greedy(n_arms =n_arms, epsilon=params.bandit$epsilon, reward_func=params.bandit$reward_func)
  #### Wishlist: add bandit algorithms besides UCB/RC/egreedy [must also update in "bandit.R"]
  } else {
    sprintf("Error: cannot call bandit (Bandit type unknown, not specified) %s",params.bandit$name)
  }
  
  # print(bandit)
  
  #Now run the search
  for (trial in times) {
    ps <- NULL
    for (arm_idx in seq(n_arms)) {
      ps <- c(ps, probability_arm(bandit, arm_idx))
    }
    # cat("trial: ", trial, "   ", "ps: ", ps, "\n")
    chosen_arm <- next_arm(bandit)
    #print(chosen_arm)
    
    search_stats <- pull_arm(chosen_arm, search_stats,block.size,SearchStrategy)
    reward       <- RewardFunction(search_stats[[chosen_arm]],block.size)
    #reward function of the form RewardFunction(new_st,block.size) 
    # unfound      <- UnfoundPrevalence(search_stats[[chosen_arm]])
    bug         <- BugsFound(search_stats[[chosen_arm]],block.size)
    
    # cat("  arm: ", chosen_arm, "  reward", reward, " unfound", unfound, "total bugs&houses", bug, "\n")
    
    # IMPORTANT: The order in which chiris are presented to the bandit from a block search is randomized 
    randomizer <- runif(block.size) 
    random.chiris <-data.frame(cbind(reward,randomizer))
    random.chiris <- random.chiris[order(randomizer),]
    # print(random.chiris)
    
    #Update the bandit for each time an arm was searched based on the update scheme for the bandit
    for (house.counter in block.size) {
      # print(random.chiris$reward[block])
      bandit <- update_bandit(bandit, chosen_arm, random.chiris$reward[house.counter])  #change bandit preferences
    }
    # cat("  preferences:", paste(bandit$preferences), "\n")
    
    #Update tables for benchmarking
    #Rewards reported are the average reward in the set of pulls
    rewards[trial] <- sum(reward)/block.size
    arms[trial]    <- chosen_arm
    #     ps_holding[trial,] <-ps
    #     mR[trial] <- bandit$mean_reward
    #     unfoundprev[trial] <- unfound
    blocking[trial] <-block.size
    bugs.houses[trial,] <-bug #uncommented to get total bugs
  }
  # colnames(ps_holding)=params.arm
  colnames(bugs.houses)=c("total.bugs.found","total.infested.houses.found") #uncommented to get total bugs
  
  #Benchmarking funcitons 
  #   results <- data.frame(ps_holding,bugs.houses,
  #                         T=times, ChosenArm=arms, BlockAvgReward=rewards, CumulativeAvgReward=cumsum(rewards*blocking), 
  #                         MeanReward=mR, UnfoundPrev.ChosenArm=unfoundprev, BlockSize=blocking,
  #                         static.rewards.available,static.bugs.available,static.houses.infested)
  #   results$CumulativeHousesFound<-cumsum(results$total.infested.houses.found)
  #   results$CumulativeBugsFound<-cumsum(results$total.bugs.found)
  # results<-sum(rewards*blocking)
  results <- data.frame(total.bugs.found=cumsum(bugs.houses[,"total.bugs.found"]),
                        total.houses.found=cumsum(bugs.houses[,"total.infested.houses.found"]),
                        zeit=times, ChosenArm=arms, houses.searched=block.size*times,
                        static.bugs.available,static.houses.infested)
  
  # results<-tail(results,1)
  rownames(results) <- NULL
  search_stats$results <- results
  
  #Create database of houses searched in order and update latitude. This allows for plotting of searches over time
  for (arm.used in 1:n_arms) { 
    order <- results[which(results$ChosenArm==arm.used),"houses.searched"]
    order <-rep(order,10)
    order<-sort(order)
    order <-cbind(order,search_stats[[arm.used]]$visited_coordinates)
    order$chosen.arm <-arm.used
    names(order) <-c("order","lat","lon","chosen.arm")
    if (arm.used ==1) {
      sites.visited <-order
    }
    else {
      sites.visited <-rbind(sites.visited,order)
    }
  }
  sites.visited$combined.lat<-(sites.visited$chosen.arm-1)*params$nrow+sites.visited$lat
  sites.visited <-sites.visited[order(sites.visited$order),]
  rownames(sites.visited) <- NULL
  search_stats$sites.visited <- sites.visited
  return(search_stats)
}