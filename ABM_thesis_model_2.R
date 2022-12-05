

##########################################
## Initial agents and time step parameters
##########################################

num_superindividuals <- 1500
Sr = 6666 ## number of individuals a superindividuals represents
num_time_steps <- 600
vmax <- 0.097083333 #(nmol/nmol C*hour)
km <- 510 #(nmol/L)
m0 <- 0.00022 #(nmole C/cell)
matrix_size <- 10000

## Set up data frame for state variables
state <- data.frame("S" = 500)
state$flow_rate <- 0.02083333
state$nutrient_inflow <- 500

##########################################
## Set up vector for agents
##########################################

agents <- matrix(NA, nrow = matrix_size, ncol = 8)
agents[1:num_superindividuals, 1] <- 1:num_superindividuals  #ID
agents[1:num_superindividuals, 2] <- 0.00033  ## cell size
agents[1:num_superindividuals, 3] <- 1 ## alive
agents[1:num_superindividuals, 4] <- state$S ## S
agents[1:num_superindividuals, 5] <- 0.00000079 ## internal resource concentration (nmol/cell)
agents[1:matrix_size, 6] <- 0.04166667 ## mumax (1/hour)
agents[1:matrix_size, 7] <- 0.002333 ## qnaught (nmol/nmol C)
agents[1:matrix_size, 8] <- sample(1:40, matrix_size, replace = TRUE, prob = c(2,rep(1,39)))
  
## make another matrix for colony information
## convection and death on colony level, then remove those individuals from agents matrix
## update colony matrix at the end of the timestep with updated values
## tracemem(agents) doesn't print if not copied
## keep past timesteps for colony dataframe
## assign() get() 
## pull mean value for colony size from papers (approx. how many cells in a colony)

#agents[, 6] <- runif(250000000, min = 1.0, max = 1.2) ## mumax w/ variability
#agents[, 7] <- runif(250000000, min = 0.00000030, max = 0.00000040) ## qnaught w/ variability

#############################################
## Set up diagnostic and plotting data frames
#############################################

output <- replicate(num_time_steps, NA, simplify = FALSE)
resource_time <- data.frame()
intracellular_time <- data.frame()
uptake_time <- data.frame()
uptake_per_individual <- data.frame()
Sr_time <- data.frame()
intermediate_intracellular<- data.frame() ##CMG TRIED ADDING THIS
intermediate_intracellular_concentration <- data.frame()
foo_time <- data.frame()

##########################################
## for loop for model run
##########################################

#profvis({
for (i in 1:num_time_steps) {
  which(!is.na(agents[,3]) & agents[,3] == 1) |> length() -> agents_before_death
  intracellular_before_death <- sum(agents[1:num_superindividuals, 5]*Sr)
  ## agent death
  random <- runif(num_superindividuals, 0, 100)
  agents[seq_len(num_superindividuals),][random <= 100*(0.5/24), 3] <- 0 #CMG CHANGED
  
  random[which(random<100*(0.5/24))]<-0
  random[which(random>=100*(0.5/24))]<-1
  agents[1:num_superindividuals,3] <- random
  
  
  ## filter by live agents 
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1)
  agents_after_death <- length(live_agents)
  intracellular_after_death <- sum(agents[live_agents, 5]*Sr)
  #agents_alive <- agents[!is.na(agents[,3]) & agents[,3] == 1,]
  
  ## determine uptake for live agents
  uptake_per_agent <- vmax*(state$S /(km + state$S))
  foo <- (state$S/(km + state$S))
  ## multiply uptake per individual by each individual's cell size
  resource_uptake_per_agent <- uptake_per_agent*agents[live_agents, 2]
  resource_uptake_all_ind <- resource_uptake_per_agent*Sr
  ## sum uptake per individual for extracellular mass balance
  uptake_total <- sum(resource_uptake_all_ind)
  
  ## incorporate uptake into agents
  if (state$S < uptake_total){
    uptake_S <- (state$S/(length(live_agents)*Sr))
    agents[live_agents, 5] <- agents[live_agents, 5] + uptake_S
    uptake_total = uptake_S
    state$S = 0 
  }
  else{
    agents[live_agents, 5] <- agents[live_agents, 5] + resource_uptake_per_agent
    state$S = state$S - uptake_total
  }
  ## update S by dilution
  state$S <- state$S - (state$flow_rate*state$S)
  
  intracellular_after_uptake_before_division <- sum(agents[live_agents, 5]*Sr,na.rm=TRUE)
  
  ## calculate mu (growth rate), convert intracellular resource to quota (nmol/nmol C)
  mu <- agents[live_agents, 6]*(1-(agents[live_agents, 7]/
                                     (agents[live_agents, 5]/agents[live_agents, 2])))
  ## update individual cell size of all agents
  agents[live_agents, 2] <- agents[live_agents, 2] + (mu*agents[live_agents, 2])
  ## if size calculated by division function is below minimum cell size (0.75), 
  ## replace with minimum cell size
  if (any(!is.na(agents[, 2]) & agents[, 2] < 0.00022)) { 
    warning("size too small")
  }
  
  ## division by cell size
  
  ##### need to add colony number!!!!
  last_ID <- tail(live_agents, n = 1)
  for (j in live_agents){
    if (agents[j, 2] > (2*m0)){
      agents[j, 2] = agents[j, 2]/2 # size update
      
      agents[j, 5] = agents[j, 5]/2 # intracellular resource update
      new <- c(last_ID + 1, agents[j, 2], 1, state$S, agents[j, 5], agents[j, 6],
               agents[j, 7])
      num_superindividuals <- num_superindividuals + 1
      total_superindividuals <- (which(!is.na(agents[,1])))
      agents[length(total_superindividuals) + 1,] <- new
      last_ID <- last_ID + 1
    }
  }
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1)
  agents_after_division <- length(agents[live_agents])
  
  ##update S with new inflow nutrients
  state$S <- state$S + (state$flow_rate*state$nutrient_inflow) 
  
  
  ## update S in agents matrix
  agents[, 4] <- state$S
  #agents_alive <- agents[!is.na(agents[,3]) & agents[,3] == 1,]
  out_per_time <- data.frame("ID" = agents[live_agents, 1], 
                             "vmax" = vmax, 
                             "km" = km,
                             "S" = state$S, 
                             "mumax" = agents[live_agents, 6],
                             "qnaught" = agents[live_agents, 7],
                             "time" = rep(i, nrow(agents[live_agents,])), 
                             "intracellular_resource" = agents[live_agents, 5],
                             "mu" = mu,
                             "cell_size" = agents[live_agents, 2],
                             "alive" = agents[live_agents, 3])
  ## update number of superindividuals and reframe agents vector for next timestep
  
  ## the adjustment recommended leaves around 50k individuals unaccounted for
  ## can get rid of this
  num_superindividuals <- length(agents[live_agents])
  if ((num_superindividuals > 2000) | (num_superindividuals < 1000)){
    Sr = round(((Sr*num_superindividuals)/1500), digits = 0)
    if (Sr <= 0){
      Sr = 1
    }
    num_superindividuals = 1500
  }
  ## filter by live agents again after division
  live_agents <- which(!is.na(agents[,3]) & agents[,3] == 1)
  intracellular_resource <- agents[live_agents, 5]
  cell_size <- agents[live_agents, 2]
  agents <- matrix(NA, nrow = matrix_size, ncol = 7)
  agents[1:num_superindividuals, 1] <- 1:num_superindividuals  #ID
  agents[1:num_superindividuals, 2] <- cell_size[1]  ## cell size
  agents[1:num_superindividuals, 3] <- 1 ## alive
  agents[1:num_superindividuals, 4] <- state$S ## S
  agents[1:num_superindividuals, 5] <- intracellular_resource[1] ## intracellular resource (nmol/cell)
  agents[1:matrix_size, 6] <- 0.04166667 ## mumax (1/hour)
  agents[1:matrix_size, 7] <- 0.002333 ## qnaught (nmol/nmol C)
  
  
  ## creating timestep data frames for diagnostics and plotting
  resource_per_time <- data.frame("time" = i, "concentration" = state$S)
  intracellular_per_timestep <- data.frame("time" = i, "intracellular_after_uptake_before_division" = intracellular_after_uptake_before_division)
  intermediate_intracellular <- data.frame("time" = i, "intracellular_before_death" = intracellular_before_death, 
                                           "intracellular_after_death" = intracellular_after_death)
  uptake_per_timestep <- data.frame("time" = i, "uptake_size" = uptake_total)
  uptake_per_individual_timestep <- data.frame("time" = i, "uptake_per_agent" = uptake_per_agent[1])
  Sr_per_timestep <- data.frame("time" = i, "Sr" = Sr, "num_superindividuals" = num_superindividuals,
                                "agents_before_death" = agents_before_death, "agents_after_death" = agents_after_death,
                                "agents_after_division" = agents_after_division)
  foo_timestep <- data.frame("time" = i, "S" = state$S, "foo" = foo)
  
  uptake_time <- rbind(uptake_time, uptake_per_timestep)
  resource_time <- rbind(resource_time, resource_per_time)
  intracellular_time <- rbind(intracellular_time, intracellular_per_timestep)
  uptake_per_individual <- rbind(uptake_per_individual, uptake_per_individual_timestep)
  Sr_time <- rbind(Sr_time, Sr_per_timestep)
  intermediate_intracellular_concentration <- rbind(intermediate_intracellular_concentration, intermediate_intracellular)
  foo_time <- rbind(foo_time, foo_timestep)
  
  #row-append to output dataframe that stores outputs for all agents and all time steps
  output[[i]] <- out_per_time
}
#})

output_df <- data.frame()

for (i in 1:num_time_steps){
  output_interem <- as.data.frame(output[[i]])
  output_df <- rbind(output_df, output[[i]])
}

intracellular_time <- cbind(intracellular_time, resource_time[,2, drop=FALSE])
intracellular_time$mass_balance <- (intracellular_time$intracellular_after_uptake_before_division+intracellular_time$concentration)

library(ggplot2)
library(dplyr)

ggplot(Sr_time, aes(x = time, y = (num_superindividuals))) + geom_point() ##USEFUL
ggplot(Sr_time, aes(x = time, y = (num_superindividuals*Sr))) + geom_point() ##USEFUL
ggplot(resource_time, aes(x = time, y = concentration)) + geom_point()
ggplot(intracellular_time, aes(x = time, y = intracellular_per_time)) + geom_point(aes(color = "intracellular concentration")) +
  geom_point(aes(y = concentration, color = "S"))
ggplot(intracellular_time, aes(x = time, y = intracellular_after_uptake_before_division)) + geom_point()##USEFUL
ggplot(intracellular_time, aes(x = time, y = mass_balance)) + geom_point()##USEFUL

ggplot(uptake_time, aes(x = time, y = uptake_size)) + geom_point()
ggplot(Sr_time, aes(x = agents_after_death, y = agents_after_division)) + geom_point()

ggplot(foo_time, aes(x = S, y = foo)) + geom_point()

plot(output_df$mu ~ output_df$q) #DOES NOT WORK??