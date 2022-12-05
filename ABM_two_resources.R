library(profvis)

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
num_agents <- 1000
num_time_steps <- 1000
flow_rate <- 0.86
nutrient_inflow <- 1.4
vmax1 <- 0.0000014
vmax2 <- 0.0000012
km1 <- 17
km2 <- 18
## mumax <- 1.1
## qnaught <- 0.00000036
m0 <- 0.00075


##### start quota at twice qnaught
## Set up vector for agents

agents <- matrix(NA, nrow = 2000000, ncol = 7)
agents[1:num_agents, 1] <- 1:num_agents  #ID
agents[1:num_agents, 2] <- 0.000075  ## cell size
agents[1:num_agents, 3] <- 1 ## alive
agents[1:num_agents, 4] <- 0.14 ## S
agents[1:num_agents, 5]
agents[1:num_agents, 6] <- 0.0000004 ## q
agents[, 7] <- runif(2000000, min = 1.0, max = 1.2) ## mumax w/ variability
agents[, 8] <- runif(2000000, min = 0.00000030, max = 0.00000040) ## qnaught w/ variability

## Set up data frame for state variables
state <- data.frame("S1" = 0.14, "S2 = 0.13")
state$flow_rate <- 0.86
state$nutrient_inflow <- 0.75
## Set up other parameter values

output <- replicate(num_time_steps, NA, simplify = FALSE)
agents_time <- data.frame()
resource_time <- data.frame()


## make data frame for time it takes to do each time step
#profvis({
for (i in 1:num_time_steps) {
  random <- runif(num_agents, 1, 100)
  ## death rate
  agents[seq_len(num_agents),][random < 15, 3] <- 0
  ## determine quota for all agents, both dead and alive
  
  quota_per_agent <- vmax*(agents[,4]/(km + agents[,4]))
  ## add quota for each agent (dead and alive) to existing quota
  ## FIX THIS FOR 2 S
  agents[, 5] <- agents[, 5] + quota_per_agent
  ## filter by live agents
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1)
  ## replace quota values less than qnaught to be qnaught
  agents[live_agents, 5] <- ifelse(agents[live_agents, 5] < 0.00000036, 
                                   0.00000036, agents[live_agents, 5])
  ## determine total quota for time step for all live agents
  quota_per_time <- sum(quota_per_agent[live_agents])
  ## determine growth rate for all agents dead and alive
  mu <- agents[, 7]*(1-(agents[, 8]/agents[, 5]))
  ## update individual cell size of all agents
  ## agents[, 2] <- (1 + mu)*agents[, 2] **not sure why there is a (1 + mu)
  agents[, 2] <- agents[, 2] + (mu*agents[, 2])
  ## if size calculated by division function is below minimum cell size (0.75), 
  ## replace with minimum cell size
  if (any(!is.na(agents[, 2]) & agents[, 2] < .00075)) { 
    warning("size too small")
  }
  agents[live_agents, 2] <- ifelse(agents[live_agents, 2] < 0.00075, 0.00075, agents[live_agents, 2])
  ## division by cell size
  for (j in live_agents){
    if (agents[j, 2] > (2*m0)){
      agents[j, 2] = m0 # size update
      agents[j, 5] = agents[j, 5]/2 # q update
      new <- c(agents[num_agents, 1] + 1, m0, 1, state$S, agents[j, 5], agents[j, 6],
               agents[j, 7])
      num_agents <- num_agents + 1
      agents[num_agents,] <- new
      
      mu[num_agents] <- mu[j]
    }
    
  }
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1) # contains row numbers of live agents
  ## update extracellular nutrient concentration state variable
  state$quota_per_time <- quota_per_time
  new_S <- state$S + apply(state, MARGIN = 1, FUN = function(x) extracellular(x))
  state <- replace(state, 1, new_S)
  if (state$S < 0){
    state$S <- 0
  }
  ## update S in agents matrix
  agents[, 4] <- state$S
  agents_alive <- agents[!is.na(agents[,3]) & agents[,3] == 1,]
  out_per_time <- data.frame("ID" = agents_alive[,1], 
                             "vmax" = vmax, 
                             "km" = km,
                             "S" = state$S, 
                             "mumax" = agents_alive[, 6],
                             "qnaught" = agents_alive[, 7],
                             "time" = rep(i, nrow(agents_alive)), 
                             "q" = agents_alive[, 5],
                             "mu" = mu[live_agents],
                             "cell_size" = agents_alive[, 2],
                             "alive" = agents_alive[,3])
  
  agent_num_per_time <- nrow(agents_alive)
  agents_per_time <- data.frame("time" = i, "number_agents" = agent_num_per_time)
  resource_per_time <- data.frame("time" = i, "concentration" = state$S)
  resource_time <- rbind(resource_time, resource_per_time)
  agents_time <- rbind(agents_time, agents_per_time)
  
  
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output[[i]] <- out_per_time
}
#})

output_df <- data.frame()

for (i in 1:111){
  output_interem <- as.data.frame(output[[i]])
  output_df <- rbind(output_df, output[[i]])
}

library(ggplot2)
ggplot(output_df, aes(x = time, y = mu)) + geom_point() + xlim(0, 45) + ylim(0, 1.15)
ggplot(output_df, aes(x = time, y = q)) + geom_point()
ggplot(agents_time, aes(x = time, y = log(number_agents))) + geom_point()
## ggplot(output_df, aes(x = time, y = S)) + geom_point()
ggplot(resource_time, aes(x = time, y = concentration)) + geom_point()