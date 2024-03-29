
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
  ## multiply vmax by 1 million
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
num_agents <- 200
num_time_steps <- 30
flow_rate <- 0.86
nutrient_inflow <- 1.4
vmax <- 0.0000014
km <- 17
mumax <- 1.1
qnaught <- 0.00000036
m0 <- 0.00075


##### start quota at twice qnaught
## Set up data frame for agents
## agents <- data.frame("ID" = 1:num_agents, "vmax" = 0.0000014, "km" = 17, 
                     ## "mumax" = 1.1, "qnaught" = 0.00000036, "cell_size" = 0.00075, 
                     ## "alive" = 1)
##todo: STORE IDENTICAL VARIABLES SEPARATELY
agents <- matrix(NA, nrow = 10000, ncol = 5)
agents[1:num_agents, 1] <- 1:num_agents  #ID
agents[1:num_agents, 2] <- 0.000075  ## cell size
agents[1:num_agents, 3] <- 1 ## alive
agents[1:num_agents, 4] <- 1.2 ## S
agents[1:num_agents, 5] <- 0 ## q

## Set up data frame for state variables
state <- data.frame("S" = 1.2)
## try with 14 for inflow and starting concentration
## agents <- cbind(agents, state)
#set seed here
#state$flow_rate <- runif(1, min = 0.16, max = 0.86)
state$flow_rate <- 0.86
state$nutrient_inflow <- 1.4
## Set up other parameter values

## todo: will want to make these long vector and then fill it in
#q <- rep(0, num_agents)
#size <- rep(0, num_agents)
output <- replicate(num_time_steps, NA, simplify = FALSE)
agents_time <- data.frame()


## make data frame for time it takes to do each time step
profvis({
for (i in 1:num_time_steps) {
  #print(i)
  ## mortality (**use dilution rate, in Hellweger paper is 0.5 per day for Droop paper)
  random <- runif(num_agents, 1, 100)
  agents[seq_len(num_agents),][random < 10, 3] <- 0 
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1) # contains row numbers of live agents
  # for (k in agents[1:num_agents,1]){
  #   random <- runif(1, min = 1, max = 100)
  #   if (random < 10.0){
  #     agents[k,3] = 0
  #   }
  # }
  ## filter out dead agents
  #agents_alive <- agents[live_agents,]
  ## update S from state data frame
  ## agents_alive$S <- state$S
  # q <- q[live_agents]
  ## sum uptake per time step to use for extracellular concentration
  ## quota_per_agent <- apply(agents_alive[, c(2, 3, 8)], MARGIN = 1, FUN = function(x) uptake(x))
  quota_per_agent <- vmax*(agents[,4]/(km + agents[,4]))
  quota_per_time <- sum(quota_per_agent[live_agents])
  #q <- q + apply(agents_alive[, c(2, 3, 8)], MARGIN = 1, FUN = function(x) uptake(x))
  ## tmp <- cbind(agents_alive, q)
  mu <- mumax*(1-(qnaught/agents[, 5]))
  ## mu <- apply(tmp[, c(4, 5, 9)], MARGIN = 1, FUN = function(x) growth(x))
  ## tmp <- cbind(tmp, mu)
  #size <- size[live_agents]
  ## size <- size + apply(tmp[, c(6, 10)], MARGIN = 1, FUN = function(x) division(x))
  agents[, 2] <-  (1 + mu)*agents[, 2]
  ## if size calculated by division function is below minimum cell size (0.75), 
  ## replace with minimum cell size
  if (any(!is.na(agents[, 2]) & agents[, 2] < .00075)) { 
    warning("size too small")
  }
  agents[live_agents, 2] <- ifelse(agents[live_agents, 2] < 0.00075, 0.00075, agents[live_agents, 2])
#  index = 1
 # browser()
  ## division by cell size
  for (j in live_agents){
    if (agents[j, 2] > (2*m0)){
      agents[j, 2] = m0 # size update
      agents[j, 5] = agents[j, 5]/2 # q update
      new <- c(agents[num_agents, 1] + 1, m0, 1, state$S, agents[j, 5])
      num_agents <- num_agents + 1
      agents[num_agents,] <- new
      #q <- append(q, q[j])
      mu[num_agents] <- mu[j]
      #size <- append(size, size[j])
      #index = index + 1
      ## write a function that can be applied row wise
      ## try putting bind further up
    }
    #index = index + 1
  }
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1) # contains row numbers of live agents
  #agents_alive <- agents[!is.na(agents[,3]) & agents[,3] == 1,]
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
  # #dataframe of outputs for each time step
  # print(agents_alive[,1])
  # print(rep(i, nrow(agents_alive)))
  # print(agents_alive)
  # print(q)
  # print(length(q))
  # print(size)
  # print(length(size))
  agents_alive <- agents[!is.na(agents[,3]) & agents[,3] == 1,]
  out_per_time <- data.frame("ID" = agents_alive[,1], 
                             "vmax" = vmax, 
                             "km" = km,
                             "S" = state$S, 
                             "mumax" = mumax,
                             "qnaught" = qnaught,
                             "time" = rep(i, nrow(agents_alive)), 
                             "q" = agents_alive[, 5],
                             "mu" = mu[live_agents],
                             "cell_size" = agents_alive[, 2],
                             "alive" = agents_alive[,3])
  
  agent_num_per_time <- nrow(agents_alive)
  agents_per_time <- data.frame("time" = i, "number_agents" = agent_num_per_time)
  agents_time<- rbind(agents_time, agents_per_time)

  
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output[[i]] <- out_per_time
}
})

## profvis package
library(ggplot2)
one <- filter(output, ID == 1)
ggplot(one, aes(x = time, y = q)) + geom_point()
ggplot(output, aes(x = time, y = mu)) + geom_point() + xlim(0, 45) + ylim(0, 1.15)
ggplot(output, aes(x = time, y = q)) + geom_point()
ggplot(agents_time, aes(x = time, y = log(number_agents))) + geom_point()
write.csv(output,"C:\\Users\\Matt\\Desktop\\Thesis\\output0.14.csv", row.names = FALSE)
read.csv("C:\\Users\\Matt\\Desktop\\Thesis\\output0.14.csv")
