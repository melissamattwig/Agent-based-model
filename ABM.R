
## Agent-based model for phytoplanton consumer resource models
## Based on Hellweger & Kianirad (2007)

library(dplyr)

###############################################
## Initial simple model for uptake and growth
###############################################

uptake <- function(x){
  vmax <- x[1]
  km <- x[2]
  s <- x[3]
  result <- vmax*(s/(km+s))
  return(result)
}

growth <- function(x){
  mumax <- x[1]
  qnaught <- x[2]
  q <- x[3]
  result <- mumax*(1-(qnaught/q))
  return(result)
}
x <- c(1.1, 0.36, 0.106)
growth(x)

num_agents <- 4
num_time_steps <- 10
agents <- data.frame("ID" = 1:num_agents, "vmax" = 1.4, "km" = 17, 
                     "s" = runif(num_agents, min = 1.3, max = 1.5),
                     "mumax" = 1.1, "qnaught" = 0.36, "Cell size" = 0.75)
q <- rep(0, num_agents)
mu <- rep(0, num_agents)
output <- data.frame()


for (i in 1:num_time_steps) { 
  q <- q + apply(agents[, c(2, 3, 4)], MARGIN = 1, FUN = function(x) uptake(x))
  tmp <- cbind(agents, q) 
  mu <- apply(tmp[, c(5, 6, 7)], MARGIN = 1, FUN = function(x) growth(x))
  #dataframe of outputs for each time step
  out_per_time <- data.frame("ID" = agents$ID, 
                             "vmax" = agents$vmax, 
                             "km" = agents$km,
                             "s" = agents$s, 
                             "mumax" = agents$mumax,
                             "qnaught" = agents$qnaught,
                             "time" = rep(i, nrow(agents)), 
                             "q" = q,
                             "mu" = mu)
  
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output <- rbind(output, out_per_time)
}


############################
### test area for division
############################
uptake <- function(x){
  vmax <- x[1]
  km <- x[2]
  s <- x[3]
  result <- vmax*(s/(km+s)) #add internal quota
  return(result)
}

growth <- function(x){
  mumax <- x[1]
  qnaught <- x[2]
  q <- x[3]
  result <- mumax*(1-(qnaught/q))
  return(result)
}

division <- function(x){
  mu <-x[1]
  m <- x[2]
  result <- mu*m
  return(result)
}

extracellular <- function(x){
  flow_rate <- x[1]
  nutrient_inflow <- x[2]
}

## set up initial conditions
num_agents <- 5
num_time_steps <- 1
agents <- data.frame("ID" = 1:num_agents, "vmax" = 0.0000014, "km" = 17, 
                     "s" = runif(num_agents, min = 1.3, max = 1.5),
                     "mumax" = 1.1, "qnaught" = 0.00000036, "cell_size" = 0.00075, 
                     "alive" = 1)
q <- rep(0, num_agents)
mu <- rep(0, num_agents)
m0 <- 0.00075
output <- data.frame()

## create data frame for global variables/state variables

for (i in 1:num_time_steps) {
  q <- q + apply(agents[, c(2, 3, 4)], MARGIN = 1, FUN = function(x) uptake(x))
  tmp <- cbind(agents, q) 
  mu <- apply(tmp[, c(5, 6, 9)], MARGIN = 1, FUN = function(x) growth(x))
  tmp <- cbind(tmp, mu)
  size <- apply(tmp[, c(7, 10)], MARGIN = 1, FUN = function(x) division(x))
  index = 1
  for (j in size){
    if (j > 0.77){
      size[index] = m0 
      new <- c(nrow(agents) + 1, 1.4, 17,
             runif(1, min = 1.3, max = 1.5), 1.1,
             0.36, m0)
      agents <- rbind(agents, new)
      q <- append(q, q[index])
      mu <- append(mu, mu[index])
      size <- append(size, size[index])
      index = index + 1
    }
  }
  #dataframe of outputs for each time step
  out_per_time <- data.frame("ID" = agents$ID, 
                             "vmax" = agents$vmax, 
                             "km" = agents$km,
                             "s" = agents$s, 
                             "mumax" = agents$mumax,
                             "qnaught" = agents$qnaught,
                             "time" = rep(i, nrow(agents)), 
                             "q" = q,
                             "mu" = mu,
                             "cell_size" = size)
                             
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output <- rbind(output, out_per_time)
}

## division rate from Hellweger & Kianirad, currently using lower value above to test model
## (2*tmp[index,7])

############################
### test area for mortality
############################

uptake <- function(x){
  vmax <- x[1]
  km <- x[2]
  S <- x[3]
  result <- vmax*(S/(km+S))
  return(result)
}

growth <- function(x){
  mumax <- x[1]
  qnaught <- x[2]
  q <- x[3]
  result <- mumax*(1-(qnaught/q))
  return(result)
}

division <- function(x){
  mu <-x[1]
  m <- x[2]
  result <- mu*m
  return(result)
}

extracellular <- function(x){
  S <- x[1]
  flow_rate <- x[2]
  nutrient_inflow <- x[3]
  quota_per_time <- x[4]
  result <- ((flow_rate*nutrient_inflow)-(flow_rate*S)-quota_per_time) # only change in quota for each time step
  return(result)
}


## Initial agents and time step parameters
num_agents <- 100
num_time_steps <- 800


##### start quota at twice qnaught
## Set up data frame for agents
agents <- data.frame("ID" = 1:num_agents, "vmax" = 0.0000014, "km" = 17, 
                     "mumax" = 1.1, "qnaught" = 0.00000036, "cell_size" = 0.00075, 
                     "alive" = 1)
## Set up data frame for state variables
state <- data.frame("S" = 1.4)
## try with 14 for inflow and starting concentration
agents <- cbind(agents, state)
#set seed here
#state$flow_rate <- runif(1, min = 0.16, max = 0.86)
state$flow_rate <- 0.86
state$nutrient_inflow <- 1.4
## Set up other parameter values
q <- rep(0, num_agents)
size <- rep(0, num_agents)
m0 <- 0.00075
output <- data.frame()
agents_time <- data.frame()

### set up dataframes to keep track of times
mortality_time <- data.frame()
functions_time <- data.frame()
division_time <- data.frame()

startTime <- Sys.time()
## make data frame for time it takes to do each time step
for (i in 1:num_time_steps) {
  ## mortality (**use dilution rate, in Hellweger paper is 0.5 per day for Droop paper)
  mortality_start <- Sys.time()
  for (k in agents$ID){
    random <- runif(1, min = 1, max = 100)
    if (random < 2.0){
      agents[k,7] = 0
    }
  }
  mortality_end <- Sys.time()
 # for (k in agents$ID){
    
 # }
  ## filter out dead agents
  agents_alive <- filter(agents, alive == 1)
  ## update S from state data frame
  agents_alive$S <- state$S
  q <- q[1:nrow(agents_alive)]
  ## sum uptake per time step to use for extracellular concentration
  functions_start <- Sys.time()
  quota_per_agent <- apply(agents_alive[, c(2, 3, 8)], MARGIN = 1, FUN = function(x) uptake(x))
  quota_per_time <- sum(quota_per_agent)
  q <- q + apply(agents_alive[, c(2, 3, 8)], MARGIN = 1, FUN = function(x) uptake(x))
  q <- ifelse(q < 0, 0, q)
  tmp <- cbind(agents_alive, q)
  mu <- apply(tmp[, c(4, 5, 9)], MARGIN = 1, FUN = function(x) growth(x))
  tmp <- cbind(tmp, mu)
  size <- size[1:nrow(agents_alive)]
  size <- size + apply(tmp[, c(6, 10)], MARGIN = 1, FUN = function(x) division(x))
  functions_end <- Sys.time()
  ## if size calculated by division function is below minimum cell size (0.75), 
  ## replace with minimum cell size
  size <- ifelse(size < 0.00075, 0.00075, size)
  index = 1
  ## division by cell size
  
  ####
  #TIME HOG
  ####
  division_start <- Sys.time()
  for (j in size){
    if (j > (2*tmp[index,6])){
      size[index] = m0
      q[index] = (q[index]/2)
      last <- tail(agents$ID, n = 1)
      new <- c(last + 1, 0.0000014, 17,
               1.1, 0.00000036, m0, 1, state$S)
      agents <- rbind(agents, new)
      agents_alive <- rbind(agents_alive, new)
      q <- append(q, q[index])
      mu <- append(mu, mu[index])
      size <- append(size, size[index])
      index = index + 1
      ## write a function that can be applied row wise
      ## try putting bind further up
    }
    index = index + 1
  }
  division_end <- Sys.time()
  ## update extracellular nutrient concentration state variable
  state$quota_per_time <- quota_per_time
  
  ############## ONLY QUOTA UPTAKE FROM TIME STEP SHOULD BE SUBTRACTED
  ## extracellular function should be sum of uptake rate across all individuals 
  new_S <- state$S + apply(state, MARGIN = 1, FUN = function(x) extracellular(x))
  state <- replace(state, 1, new_S)
  if (state$S < 0){
    state$S <- 0
  }
  #state <- ifelse(state$S < 0, 0, state)
  #dataframe of outputs for each time step
  out_per_time <- data.frame("ID" = agents_alive$ID, 
                             "vmax" = agents_alive$vmax, 
                             "km" = agents_alive$km,
                             "S" = agents_alive$S, 
                             "mumax" = agents_alive$mumax,
                             "qnaught" = agents_alive$qnaught,
                             "time" = rep(i, nrow(agents_alive)), 
                             "q" = q,
                             "mu" = mu,
                             "cell_size" = size,
                             "alive" = agents_alive$alive)
  
  agent_num_per_time <- nrow(agents_alive)
  agents_per_time <- data.frame("time" = i, "number_agents" = agent_num_per_time)
  agents_time<- rbind(agents_time, agents_per_time)
  
  mortality_per_time <- data.frame('time' = i, "start" = mortality_start, "end" = mortality_end)
  mortality_time <- rbind(mortality_time, mortality_per_time)
  
  functions_per_time <- data.frame('time' = i, "start" = functions_start, "end" = functions_end)
  functions_time <- rbind(functions_time, functions_per_time)
  
  division_per_time <- data.frame('time' = i, "start" = division_start, "end" = division_end)
  division_time <- rbind(division_time, division_per_time)
  
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output <- rbind(output, out_per_time)
}
endTime <- Sys.time()

library(ggplot2)
one <- filter(output, ID == 1)
ggplot(one, aes(x = time, y = q)) + geom_point()
ggplot(output, aes(x = time, y = mu)) + geom_point() + xlim(0, 45) + ylim(0, 1.15)
ggplot(output, aes(x = time, y = q)) + geom_point()
ggplot(agents_time, aes(x = time, y = number_agents)) + geom_point()
write.csv(output,"C:\\Users\\Matt\\Desktop\\Thesis\\output0.14.csv", row.names = FALSE)
read.csv("C:\\Users\\Matt\\Desktop\\Thesis\\output0.14.csv")
