##' Lower level helper functions NOT specific to any search algo
##' 

#TODO: set working directory to point to simple_bandit repository or create an environment variable SIMPLE_BANDIT

if(Sys.getenv("SIMPLE_BANDIT") != "") {
  setwd(Sys.getenv("SIMPLE_BANDIT"))  
} else if(Sys.getenv("USER") == "....") {
  setwd("~/simple_bandit")
}

#SET OUTPUT PATHS
if(! file.exists("output")) {
  dir.create("output")  
}

#Simple functions
#Log10 
log10plus <- function(x) {
  log10(1+x)
}

#count positive
countpositive <- function(x) {
  length(which(x>0))
}

#timing
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  toc - tic
}

TimeNow <- function() {
  x <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")
  return(x)
}


time_stamp <- function(prefix="", suffix="", outmsg=TRUE) {
  #generates a unique time stamp for every data generated
  t <- format(Sys.time(), "%Y-%m-%d__%H-%M-%S");
  s <- as.integer(runif(1, max=1000))
  filename <- paste(prefix, t, s, suffix, sep="")
  if (outmsg) {
    print(filename)
  }
  return(filename)
}

##' INFESTATION RELATED FUNCTIONS

##' Infestation Generator
##' matrix with 0=not infested, x>0= infested with x bugs
##' 1. some houses are initially infested.  
##' 2. each turn there is a probability of starting a new infestation
##  assumptions: growth is linear in the infestation size (future: use Hussl and Riker growth models)
##' @params = parameters of the infestation (as list)
    'must include: 
    ncol: number of columns for the matrix
    nrow: number of rows for the matrix
    num_starts: how many places you want the infestation to start
    p_hop: probability of hop 
    p_skip: probability of skip
    p_jump: probability of jump
    pMovePerDay = moving per day 
    pClearancePerDay = spontaneous disappearance per day (usually 0)
    '
##' @params = prevalence_rule: what the ulimate prevalence of the infestation should be
##' This will build an infestation until it reaches those parameters
GenerateInfestation <- function(params, prevalence_rule) {
  print("generating infestation ...")
  num_starts <- params$num_starts
  if(num_starts > prevalence_rule*params$nrow*params$ncol) {
    print(sprintf("Warning: num starts (%.0f) exceeds the target prevalence (%.3f) - reducing.", num_starts, prevalence_rule))
    num_starts <- as.integer(prevalence_rule*params$nrow*params$ncol)
  }
  MakeMove <- function(name, house, nrow, ncol) {
    if(name == "hop") {
      stepsize <- 1
    } else if(name == "skip") {
      stepsize <- 2
    } else {
      stopifnot(name == "jump")
      stepsize <- 3
    }
    new_house <- house + round(runif(2)*stepsize)*sign(runif(2)-c(0.5, 0.5))
    ##wraparound - toroidal grid
    #new_house[1] <- 1 + ((new_house[1]-1)%%nrow)
    #new_house[2] <- 1 + ((new_house[2]-1)%%ncol)
    new_house$lat <- min(nrow, max(1, new_house$lat))
    new_house$lon <- min(ncol, max(1, new_house$lon))
    return(new_house)
  }
  
  #CREATE INITIAL INFESTATION
  infested_coords <- data.frame()
  if (num_starts >0) {
    for(trial in 1:params$num_starts) {
      new_house <- c(sample(1:params$nrow, 1), sample(1:params$ncol, 1))
      infested_coords <- rbind(infested_coords, new_house)
    }
    names(infested_coords) <- c("lat", "lon")
  }
  
  infested_coords_prime <- infested_coords
  
  #create prevalence-based infestation: 
  #This infestation stops once a certain prevalence is reached
  prevalence=0 #intialize prevalence 
  while(prevalence<prevalence_rule)  {
    if (num_starts==0) {
      break
    }
    #clearing
    num_infested <- dim(infested_coords)[1]
    num_cleared <- rbinom(1, num_infested, params$pClearancePerDay)
    infested_coords_prime <- infested_coords[sample.int(num_infested, num_infested-num_cleared),]
    
    #could even add breeding ...
    #now move
    num_infested <- dim(infested_coords_prime)[1]
    num_spreading <- rbinom(1, num_infested, params$pMovePerDay)
    for(house_num in sample.int(num_infested, num_spreading)) {
      house <- infested_coords_prime[house_num,]
      r <- runif(1)
      if (r<params$p_hop) {
        new_house <- MakeMove("hop", house, params$nrow, params$ncol)
      } else if (r < params$p_hop + params$p_skip) {
        new_house <- MakeMove("skip", house, params$nrow, params$ncol)
      } else {
        new_house <- MakeMove("jump", house, params$nrow, params$ncol)
      }
      infested_coords_prime <- rbind(infested_coords_prime, new_house)
    }
    infested_coords <- infested_coords_prime
    
    #CALCULATE PREVALENCE (#INFECTED HOUSES/TOTAL HOUSES) to update
    prevalence <- ((dim (unique(infested_coords))[1])/(params$nrow*params$ncol))
    #print(prevalence)
  }
  
  #print(infested_coords)
  #Count the number of infested houses
  num_infested <- dim(infested_coords)[1]
  m <- matrix(0, nrow=params$nrow, ncol=params$ncol)
  if (num_starts >0) {
    for(house_num in 1:num_infested) {
      house <- unlist(infested_coords[house_num,])
      #print(house)
      #weighted
      m[house[1], house[2]] <- m[house[1], house[2]] + 1
      
      #0,1
      #m[house[1], house[2]] <- 1
    }
  }
  return(m)
}

#Generate Infestation Stats:
InfestationStats <- function(infestations=NULL) {
  total.houses.infested   <- as.double(Reduce("+",lapply(infestations,countpositive))) #total infested for benchmarking
  total.bugs.available    <- Reduce("+",lapply(infestations,sum)) #total bugs on chart
  total.rewards.available <- Reduce("+",lapply(lapply(infestations,log10plus),sum)) #total rewards if reward=log10(1+bugs)
  infest.stats <- cbind(total.houses.infested,total.bugs.available,total.rewards.available)
  colnames(infest.stats)<- c("total.houses.infested","total.bugs.available","total.rewards.available")
  colnames(infest.stats)
  return(infest.stats)
}


##' Output an image showing an infestation (assumed grid)
##'   optionally, use color to show which house was visited 
ViewInfestation <- function(infestion, visitSequence=NULL) {
  image(infestion, col=gray.colors(30, start=0.95, end=0.05))
  if(! is.null(visitSequence)) {
    #uses the color to indicate the sequence of visits
    browser()
  }
  fpath <- paste0("output/infestation", TimeNow(), ".png")  
  #quartz.save(fpath)
  ggsave(fpath)
}

#Additional method for viewing an infestation (that matches infestation matrix)
ViewInfestationGrid <-function (infestation) {
  data <- data.frame(which(infestation != 0,arr.ind=TRUE))
  data$bugs <-infestation[cbind(data$row,data$col)]
  xmax <-dim(infestation)[2]
  ymax <- dim(infestation)[1]
  pl <- ggplot(data=data,aes(x=col,y=row,col=bugs)) +geom_point(size=1.5) + coord_fixed(ratio=.5)
  pl <- pl +labs(color="# Bugs")  
  pl <- pl +scale_color_continuous (low = "yellow", high = "red", space = "Lab", na.value = "grey50")
  pl <- pl + xlim(0,xmax) +ylim(0,ymax)
  pl <- pl + labs(x="Longitude", y="Latitude", title="Infestation")
  pl <-pl +theme(legend.position = "left")
  print(pl)
  return(pl)
}




