#####################################################################

#Implementation of Bandit Heurestics

#1. Reinforcement comparison methods Sutton and Barto (1998) 
#2. UCB family of algorithms by Auer, Cesa-Bianchi & Fisher (2002)

# version 2017-06-27 by Sasha Gutfraind, Sarah Nutman
#   based on some testing code from White et al.
# ref: "Algorithms for the multi-armed bandit problem" / Precup, Kuleshov

##################
#  USAGE GUIDE
##################

# Each bandit needs:
#   -Initialize bandit function: initialize_[BANDIT_NAME]: initializes
#   -Next arm function: next_arm_[BANDIT_NAME] : pulls the bandit
#   -Probability arm function probability_arm_[BANDIT NAME]: gives the probabilities of pulling each arm
#   -Update bandit: update_bandit_[BANDIT_NAME]: function to update the bandit based on rewards
# 
# Then: 
#   Bandit must be added into generalized Next_arm, probability_arm, update_bandit function
#   (set generalized functions to call the named bandit functions)
#     
# Benchmarking function:
#   gaussian_reward_toy: determines rewards

###########################
# Example of usage:
#bandit <- initialize_[BANDIT_NAME](<args>) #must be called before any other function. Depending on type of bandit, 
#     certain additional parameters might be required (e.g. learning rate, etc)
#cycle #1
#arm_idx <- next_arm(bandit)
#r <- gaussian_reward_toy(my_system,arm_idx)
#bandit <- update_bandit(bandit,arm_idx,r)
#cycle #2
#arm_idx <- next_arm(bandit)
#r <- gaussian_reward_toy(system,arm_idx)
#bandit <- update_bandit(bandit,arm_idx,r)
#...etc

##################
#  note
##################
# arm is the same as strategy.  arms are numbered 1..n_arms

### NOTE: 
### If you add bandits to this package, you must also update the function BanditSearchGridBlocks in "HighLevelSearchAlgos.R"

#####################################################################

library(reshape2)
library(plyr)
library(ggplot2)

timeNow <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")

#RC BANDIT
#Reinforcement comparison bandit
#learning_rate=in [0,1] increases the tendency to prefer recently successful strategies
#discount_factor in [0,1];  greater for more restless bandits
# Bandit is updated with
# F=preferences
# r=reward
# R=mean reward
# A=discount factor
# B=learning rate
# 
# F(t+1)=F(t)+B(r(t)-R(t)) <- update first
# R(t+1)=A(r(t))+(1-A)(R(t))
# '

initialize_rc <- function (n_arms, learning_rate= NULL, discount_factor=NULL, preferences=NULL) {
  bandit <- list("name"="rc", "n_arms"=n_arms, learning_rate=learning_rate, discount_factor=discount_factor)
  
  if (sum(preferences) == 0) {
    bandit$preferences <- rep(0, n_arms)
  } else {
    bandit$preferences <-preferences
  }
  bandit$mean_reward <- 0.
  #print(bandit)
  return(bandit)
}

next_arm_rc <- function(bandit) {
  #returns the randomly-drawn next arm
  p_vals <- exp(bandit$preferences)
  return(sample.int(bandit$n_arms, size=1, replace=TRUE, prob=p_vals))
}


probability_arm_rc <- function(bandit, arm_idx) {
  #compute the probability for an arm
  return(exp(bandit$preferences[arm_idx]) / sum(exp(bandit$preferences)))
}

update_bandit_rc <- function(bandit, chosen_arm, reward) {
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\t", bandit$preferences, "\n")
  #cat("ZZZ  preference: ", bandit$preferences[chosen_arm], "  arm: ", chosen_arm, "   reward: ", reward, "  mean_reward: ", bandit$mean_reward, "\n")
  bandit$preferences[chosen_arm] <- bandit$preferences[chosen_arm]                + bandit$learning_rate   * (reward - bandit$mean_reward)
  bandit$mean_reward             <- (1-bandit$discount_factor)*bandit$mean_reward + bandit$discount_factor * reward 
  #bandit$mean_reward <- bandit$mean_reward + bandit$discount_factor*(reward - bandit$mean_reward)
  # cat("ZZZ  preference: ", bandit$preferences[chosen_arm], "\n")
  # cat("ZZZ  preference: ", bandit$preferences, "\n")
  # sum <- sum(exp(bandit$preferences))
  # cat("ZZZ  prob: ", exp(bandit$preferences)/sum, "\n")
  return(bandit)
}


#Upper Confidence Bound (1) bandit
#UCB
initialize_ucb1 <- function(n_arms, reward_func) {
  bandit <- list("name"="ucb1", "n_arms"=n_arms)
  bandit$mean_reward <- rep(0, n_arms)
  
  if (is.null(reward_func)) {
    for(arm_idx in seq(n_arms)) {
      bandit$mean_reward[arm_idx] <- 0
    }
  }
  else {
    for(arm_idx in seq(n_arms)) {
      bandit$mean_reward[arm_idx] <- reward_func(arm_idx, n_arms)
      }
    }
  bandit$num_pulls   <- rep(1,n_arms)
  #print(bandit)
  return(bandit)
}

next_arm_ucb1 <- function(bandit) {
  #returns the next arm.  this is a deterministic algorithm
  total_pulls <- sum(bandit$num_pulls)
  values <- list()
  for(arm_idx in seq(bandit$n_arms)) {
    values <- c(values, bandit$mean_reward[arm_idx] + sqrt( (2*log(total_pulls))/bandit$num_pulls[arm_idx] ))  
  }
  best_arm <- which.max(values)  #in cases of ties, it returns the first arm
  #print(paste(values))
  #cat("best arm: ", best_arm, "\n")
  return(best_arm)
}

probability_arm_ucb1 <- function(bandit, arm_idx) {
  #compute the probability for an arm
  cat("WARNING: this is a deterministic algorithm\n")
  best_arm <- next_arm_ucb1(bandit)
  if(arm_idx == best_arm){
    return(1.0)
  } else {
    return (0.0)
  }
}

update_bandit_ucb1 <- function(bandit, chosen_arm, reward) {
  #cat("ZZZ:", paste(bandit$num_pulls), "\n")
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\n")
  #TODO: should the first term be multiplied by (1-bandit$learning_rate)
  np <- bandit$num_pulls[chosen_arm]
  bandit$mean_reward[chosen_arm] <- (np*bandit$mean_reward[chosen_arm] + reward)/(np+1)
  bandit$num_pulls[chosen_arm]   <- np + 1
  #cat("ZZZ:", paste(bandit$num_pulls), "\n")
  #cat("ZZZ:", chosen_arm, reward, bandit$mean_reward, "\n")
  return(bandit)
}

#the test reward function: the reward is stochastic, often 0, and when positive averages
#U(arm1) = 1
#U(arm2) = 2
#...
#true_best_arm = n_arms
gaussian_reward_toy <- function(arm_idx, n_arms) { 
  stopifnot(arm_idx >= 1)
  stopifnot(arm_idx <= n_arms)
  if(runif(1) < 0.5) {
    return(0)
  } else {
    return(rnorm(1, mean=arm_idx));
  }
  return(rnorm(1, mean=arm_idx));
}


#General Bandit Functions

next_arm <- function(bandit) {
#returns the randomly-drawn next arm
  if(bandit$name=="rc") {
    return(next_arm_rc(bandit));
  } else {if (bandit$name=="ucb1") {
    return(next_arm_ucb1(bandit));
  } # Add in additional bandits here
    else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } }
}  


probability_arm <- function(bandit, arm_idx) {
#returns the current probability of selecting arm with index arm_idx
  if(bandit$name=="rc") {
    return(probability_arm_rc(bandit, arm_idx));
  } else {if (bandit$name=="ucb1") {
    return(probability_arm_ucb1(bandit, arm_idx));
  } # Add in additional bandits here
    else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } } 
}

update_bandit <- function(bandit, chosen_arm, reward) {
#called after the reward was returned
  if(bandit$name=="rc") {
    return(update_bandit_rc(bandit, chosen_arm, reward))
  } else {if (bandit$name=="ucb1") {
    return(update_bandit_ucb1(bandit, chosen_arm, reward))
  } # Add in additional bandits here
    else {
    sprintf("Error: cannot call bandit: %s",bandit$name);
  } } 
}  

