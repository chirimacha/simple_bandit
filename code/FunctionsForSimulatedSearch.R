'
FUNCTIONS FOR COMPARING RING AND BANDIT ON SIMULATED DATA
-- creates and subdivides infestations and compares running the bandit vs. doing global ring search
   Bandit search: bandit selects an arm and ring search (battleship) conducted in that arm with results reported 
   to the bandit (and arm pulled) every X times 
   Global search: Conduct ring search on a infestation created by combining the arms of the bandit into one large infestation
-- Calls "highLevelSearchAlgos.R" which contains:
    1) call to the infestation generator (create simulated data)
    2) the code to conduct a ring search on simulated data
    3) call to bandit parameters
'
#https://cran.r-project.org/web/packages/neldermead/neldermead.pdf
#library(doParallel)
library(doRNG)
library(stringr)
library(ggplot2)
library(pROC)
library(sm)
library(reshape2)
library(snowfall)

#TO DO: create an environmental variable that points to your simple_bandit repository
setwd(Sys.getenv("SIMPLE_BANDIT")) 

#Call bandit simulations
source("code/helperFunctions.R")
source("code/highLevelSearchAlgos.R")

###' ZCompareGlobalVsBandit = 
###' run a comparison of Global vs. Bandit results
###' 
###' @params.bandit = parameters of the bandit, including type of bandit. For RC: c(learning,discount) 
###' @bandit.time = number of times the bandit will search
###' @bandit.block.size= number of searches to be done in an area before returning results to the bandit
###' NOTE: bandit.time*bandit.block.size= max.global 
###' @params.grid = infestation and search parameters [see below for specific parameters to be set]
###' @params.arm = parameters for each arm of the bandit [these arms are then combined to make the global search]
###' @num.replications = number of times to repeat the experiment (i.e. compare bandit and global search)
###' @cpu.cores.needed = 1 for serial, and > 1 for parallel
ZCompareGlobalVsBandit <- function(params.bandit, bandit.time, bandit.block.size, 
                                                params.grid,  params.arm, num.replications,
                                                SearchStrategy,
                                                altAlgo, cpu.cores.needed=1) {
  compareAllReplications <- c()
  if(cpu.cores.needed==1) {
    for (r in seq(num.replications)) {
        cat("replication", r, "\n")
        compare <- ZCompareGlobalVsBandit_Helper(params.bandit=params.bandit, 
                                                 bandit.time=bandit.time, 
                                                 bandit.block.size=bandit.block.size, 
                                                 params.grid=params.grid,  
                                                 params.arm=params.arm,
                                                 SearchStrategy=SearchStrategy,
                                                 altAlgo = altAlgo) 
        
        compareAllReplications <- rbind(compareAllReplications, compare)
    }
  } else {
    print("Running in parallel...")
    sfInit(parallel=TRUE, cpus=cpu.cores.needed)
    #sfLibrary("reshape2")
    sfSource("code/helperFunctions.R")
    sfSource("code/highLevelSearchAlgos.R")    
    #sfExportAll(debug=T)
    sfExport("ZCompareGlobalVsBandit_Helper")  # #and all its arguments..)
    sfExport(list=names(formals(fun = sys.function(sys.parent()))))
    
    returnObjectList <-sfLapply(seq(num.replications), function(x){
      ZCompareGlobalVsBandit_Helper(params.bandit=params.bandit, 
                                    bandit.time=bandit.time, 
                                    bandit.block.size=bandit.block.size, 
                                    params.grid=params.grid,  
                                    params.arm=params.arm,
                                    SearchStrategy=SearchStrategy,
                                    altAlgo = altAlgo                                                                                    
    )})    
    sfStop()
    compareAllReplications <- do.call(rbind,returnObjectList)
  }
  compare.fin <-list(
    compare = compareAllReplications,
    best.arm = names(params.arm)[which.max(params.arm)]
  )
  print(t(compareAllReplications))
  return(compare.fin)
}

##' Generate single landscape and run the algorithm
ZCompareGlobalVsBandit_Helper <- function(params.bandit, bandit.time, bandit.block.size, 
                                   params.grid,  params.arm,
                                   SearchStrategy,
                                   altAlgo) {
  infestations<- lapply(params.arm, GenerateInfestation, params=params.grid)  
  stopifnot(bandit.time * bandit.block.size <= prod(dim(infestations[[1]]))*length(infestations)/2)
  
  print("A: Bandit search..")
  banditBlockSearch <- BanditSearchGridBlocks(test.time=bandit.time, params=params.grid,
                                              infestations=infestations,
                                              block.size=bandit.block.size,
                                              params.bandit=params.bandit,
                                              SearchStrategy=SearchStrategy)
  banditBlockSearch.stats <- tail(banditBlockSearch$results,1)

  print("B: alternative algorithm search..")
  #typically global = search based on global landscape with ring search around any detected infestation
  infestations.combined <- do.call(rbind,infestations) #combine infestations (i.e. arms) into one big area
  max_actions  <- bandit.time*bandit.block.size
  globalSearch <- altAlgo(infestations.combined,st=NULL,max_actions=max_actions, params=params.grid) 
  globalSearch.stats <- tail(globalSearch$running_stats, 1)[c("total_bugs","total_found","total_visited")]
  
  #show how frequently sites in each area were visited
  globalSearch.grid <- do.call(rbind, lapply(1:length(infestations), function(i){infestations[[i]]*0 + i}))  #insert i in the area of arm i
  globalSearch.arms.visited <- as.vector(globalSearch.grid * (globalSearch$visited > 0)*1)
  globalSearch.arms.visited <- table(globalSearch.arms.visited[which(globalSearch.arms.visited>0)])
  
  compare <- data.frame(globalSearch.bugsfound   = globalSearch.stats$total_bugs, 
                        globalSearch.housesfound = globalSearch.stats$total_found,
                        globalSearch.housesvisited=sum((globalSearch$visited > 0)),
                        banditBlockSearch.bugsfound   = banditBlockSearch.stats$total.bugs.found, 
                        banditBlockSearch.housesfound = banditBlockSearch.stats$total.houses.found,
                        GENERAL.bugsavailable   = banditBlockSearch.stats$static.bugs.available, 
                        GENERAL.housesinfested  = banditBlockSearch.stats$static.houses.infested)
  for(arm_idx in 1:length(infestations)) {
    armName <- names(params.arm)[arm_idx]
    compare[paste0("visitsArm_", armName, "_globalSearch")]      <- globalSearch.arms.visited[arm_idx]  
    compare[paste0("visitsArm_", armName, "_banditBlockSearch")] <- sum(banditBlockSearch[[armName]]$visited>0)
  }
  row.names(compare) <- NULL
  return(compare)
}


#create graph showing the sensitivity of the bandit and the global search
# sensitivity, i.e. ability to detect most of the infested houses (total infested houses found/total infested houses)
ZGraphResults <- function(data, params.grid, params.arm) {
  houses <- params.grid$nrow * params.grid$ncol * length(params.arm)
  banditBlockSearch.sensitivity <- (data$compare$banditBlockSearch.housesfound/data$compare$GENERAL.housesinfested)
  globalSearch.sensitivity      <- (data$compare$globalSearch.housesfound     /data$compare$GENERAL.housesinfested)
  #housesavailable == total infested houses in the landscape
  
  wtest <- tryCatch(wilcox.test(banditBlockSearch.sensitivity, globalSearch.sensitivity),
                    error=function(e){print(e); return("")})  #Wilcoxon signed-rank test on sensitivities
  print(wtest)
  
  sensitivity <-data.frame(bandit=banditBlockSearch.sensitivity, global=globalSearch.sensitivity)
  sensitivity <- rbind(data.frame(sensitivity=banditBlockSearch.sensitivity, algo="bandit"),
                       data.frame(sensitivity=globalSearch.sensitivity,      algo="global"))

  totalSearches <- mean(data$compare$globalSearch.housesvisited)
  plotInfo <- paste("Houses:",houses,"Houses Searched:",totalSearches,"\n",
                    "Arms:",str_c(params.arm,collapse=','), sep=" ")
  pl <- ggplot() 
  pl <- pl+ geom_point(data=sensitivity, aes(x=algo, y=sensitivity, color=algo))
  pl <- pl+ ylim(c(0,1)) + xlab("Search") + ylab("Sensitivity")  
  pl <- pl +ggtitle(plotInfo)
  print(pl)
  return(pl)
}


#Optional: Improved formatting for result graphing
ZGraphResults.Formated <- function(data, params.grid, params.arm, axis=c(0,1),info="OFF") {
  houses <- params.grid$nrow * params.grid$ncol * length(params.arm)
  banditBlockSearch.sensitivity <- (data$compare$banditBlockSearch.housesfound/data$compare$GENERAL.housesinfested)
  globalSearch.sensitivity      <- (data$compare$globalSearch.housesfound /data$compare$GENERAL.housesinfested)
  #housesavailable == total infested houses in the landscape
  
  wtest <- tryCatch(wilcox.test(banditBlockSearch.sensitivity, globalSearch.sensitivity),
                    error=function(e){print(e); return("")})
  print(wtest)
  
  sensitivity <-data.frame(bandit=banditBlockSearch.sensitivity, global=globalSearch.sensitivity)
  sensitivity <- rbind(data.frame(sensitivity=banditBlockSearch.sensitivity, algo="Bandit"),
                       data.frame(sensitivity=globalSearch.sensitivity,      algo="Global"))
  
  totalSearches <- mean(data$compare$globalSearch.housesvisited)
  
  plotInfo <- paste("Houses:",houses,"Houses Searched:",totalSearches,
                    "Arms:",str_c(params.arm,collapse=','), sep=" ")
  title <-"Sensitivity of Search Algorithms"
  pl <- ggplot() 
  pl <- pl+ geom_point(data=sensitivity, aes(x=algo, y=sensitivity, color=algo))
  pl <- pl+ ylim(axis) + xlab("Search Algorithm") + ylab("Sensitivity")
  # if (! is.null(info)) {
  pl <- pl +ggtitle(bquote(atop(.(title), atop(italic(.(plotInfo), "")))))
  #   }
  #   else {
  #     pl <-pl+ggtitle(bquote(.(title)))
  #   }
  # pl <- pl +ggtitle("Sensitivity of Search Algorithms",subtitle=bquote(.(plotInfo)))
  pl <- pl +labs(color = "Search Algorithm")
  
  pl <-pl +theme_bw()
  print(pl)
  return(pl)
}


##' Plot the algo for POA (percent of optimal action) and sensitivity
ZGraphResults2D <- function(data, params.grid, params.arm, searches) {
  houses <- params.grid$nrow * params.grid$ncol * length(params.arm)
  banditBlockSearch.sensitivity <- (data$compare$banditBlockSearch.housesfound/data$compare$GENERAL.housesinfested)
  globalSearch.sensitivity      <- (data$compare$globalSearch.housesfound     /data$compare$GENERAL.housesinfested)
  #housesavailable == total infested houses in the landscape
  
  optArmName <- names(which.max(params.arm))
  #algNames <- c("globalSearch", "banditBlockSearch")
  visitsData <- data$compare[,grep("visitsArm", names(data$compare))]
  visit.globalSearch      <- visitsData[,grep("globalSearch", names(visitsData))]
  visit.banditBlockSearch <- visitsData[,grep("banditBlockSearch", names(visitsData))]
  POA.globalSearch      <- visit.globalSearch[,grep(optArmName, names(visit.globalSearch))]/rowSums(visit.globalSearch)
  POA.banditBlockSearch <- visit.banditBlockSearch[,grep(optArmName, names(visit.banditBlockSearch))]/rowSums(visit.banditBlockSearch)
  
  results <- data.frame(               sensitivity=banditBlockSearch.sensitivity, POA=POA.banditBlockSearch, alg="banditBlockSearch")
  results <- rbind(results, data.frame(sensitivity=globalSearch.sensitivity,      POA=POA.globalSearch,      alg="globalSearch"))

  pl <- qplot(data=results, x=sensitivity, y=POA, color=alg, shape=alg)# + geom_point()
  pl <- pl + xlim(c(0.5, 1.0)) + ylim(c(0.5, 1.0)) + ylab("Probability of Optimal Action")
  print(pl)
  
  ggsave(time_stamp("output/global_vs_bandit_2D", ".png"), height=9, width=9)
  
  return(pl)  
}
  

#Optional functions that are used to see how the bandit changes with changing sensitivity in arms
ZScanSensitivityGlobalVsBandit <- function(params.bandit, bandit.time, bandit.block.size, 
                                           params.grid,  params.arms.static, params.arms.dynamic,
                                           num.replications,
                                           altAlgo = RingSearch) {
  
  overallResults <- c()
  armNameDynamic <- names(params.arms.dynamic)[1]
  for(prev in params.arms.dynamic[,1]) {
    params.arms <- params.arms.static
    params.arms[[armNameDynamic]] <- prev
    res <- ZCompareGlobalVsBandit(params.bandit, bandit.time, bandit.block.size, 
                           params.grid,  params.arms,
                           num.replications,
                           altAlgo = RingSearch,
                           cpu.cores.needed=4)
    banditBlockSearch.sensitivity <- (res$compare$banditBlockSearch.housesfound/res$compare$GENERAL.housesinfested)
    globalSearch.sensitivity      <- (res$compare$globalSearch.housesfound     /res$compare$GENERAL.housesinfested)
    
    overallResults <- rbind(overallResults, data.frame(vals=banditBlockSearch.sensitivity, prev=prev, algo="banditBlockSearch"))
    overallResults <- rbind(overallResults, data.frame(vals=globalSearch.sensitivity,      prev=prev, algo="globalSearch"))
    #df<-data.frame(banditBlockSearch=as.matrix(summary(banditBlockSearch)),
    #           globalSearch=as.matrix(summary(globalSearch)))
  }
  return(overallResults)
}


ZScanSensitivityGlobalVsBanditParallel <- function(params.bandit, bandit.time, bandit.block.size, 
                                           params.grid,  params.arms.static, params.arms.dynamic,
                                           num.replications, cpu.cores.needed,
                                           altAlgo = RingSearch) {
  
  overallResults <- c()
  armNameDynamic <- names(params.arms.dynamic)[1]
  sfInit(parallel=TRUE, cpus=cpu.cores.needed)
  #sfLibrary("reshape2")
  sfSource("code/helperFunctions.R")
  sfSource("code/highLevelSearchAlgos.R")    
  #sfExportAll(debug=T)
  sfExport("ZCompareGlobalVsBandit")  # #and all its arguments..)
  sfExport("ZCompareGlobalVsBandit_Helper")  # #and all its arguments..)
  sfExport(list=names(formals(fun = sys.function(sys.parent()))))
  
  returnObjectList <-sfLapply(rep(params.arms.dynamic[,1], num.replications), function(prev){
    results <- c()
    params.arms <- params.arms.static
    params.arms[[armNameDynamic]] <- prev
    res <- ZCompareGlobalVsBandit(params.bandit, bandit.time, bandit.block.size, 
                                  params.grid,  params.arms,
                                  num.replications,
                                  altAlgo = RingSearch,
                                  cpu.cores.needed=1)
    banditBlockSearch.sensitivity <- (res$compare$banditBlockSearch.housesfound/res$compare$GENERAL.housesinfested)
    globalSearch.sensitivity      <- (res$compare$globalSearch.housesfound     /res$compare$GENERAL.housesinfested)
    results <- rbind(results, data.frame(vals=banditBlockSearch.sensitivity, prev=prev, algo="banditBlockSearch"))
    results <- rbind(results, data.frame(vals=globalSearch.sensitivity,      prev=prev, algo="globalSearch"))
    return(results)
  })
  overallResults <- do.call(rbind,returnObjectList)
  return(overallResults)
}

  
