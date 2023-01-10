

##########################################
## Initial individuals and time step parameters
##########################################

num_individuals <- 45000
#Sr = 6666 ## number of individuals a individuals represents
num_time_steps <- 100
vmax <- 0.097083333 #(nmol/nmol C*hour)
km <- 510 #(nmol/L)
m0 <- 0.00022 #(nmole C/cell)
individual_matrix_size <- 100000
colony_matrix_size <- 1000

## Set up data frame for state variables
state <- data.frame("S" = 500)
state$flow_rate <- 0.02083333
state$nutrient_inflow <- 500

##########################################
## Set up colony for individuals 
##########################################

individuals <- matrix(NA, nrow = individual_matrix_size, ncol = 8)
individuals[1:num_individuals, 1] <- 1:num_individuals  #ID
individuals[1:num_individuals, 2] <- 0.00033  ## cell size
individuals[1:num_individuals, 3] <- 1 ## alive
individuals[1:num_individuals, 4] <- state$S ## S
individuals[1:num_individuals, 5] <- 0.00000079 ## internal resource concentration (nmol/cell)
individuals[1:individual_matrix_size, 6] <- 0.04166667 ## mumax (1/hour)
individuals[1:individual_matrix_size, 7] <- 0.002333 ## qnaught (nmol/nmol C)
individuals[1:individual_matrix_size, 8] <- sample(1:30, individual_matrix_size, replace = TRUE, 
                                                   prob = c(1.5, 0.75, rep(1,28))) ## colony id
  
## make another matrix for colony information
## convection and death on colony level, then remove those individuals from individuals matrix
## update colony matrix at the end of the timestep with updated values
## tracemem(individuals) doesn't print if not copied
## keep past timesteps for colony dataframe
## assign() get() 
## pull mean value for colony size from papers (approx. how many cells in a colony)

#individuals[, 6] <- runif(250000000, min = 1.0, max = 1.2) ## mumax w/ variability
#individuals[, 7] <- runif(250000000, min = 0.00000030, max = 0.00000040) ## qnaught w/ variability

##########################################
## Set up matrix for colony
##########################################

colony <- data.frame()

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
## functions created for model
##########################################

numColonies <- function(individuals){
  num_colonies <- length(unique(individuals[, 8]))
  iterator = 1
  colony_id <- c()
  colony_population <- c()
  colony_live <- c()
  for (i in 1:num_colonies){
    while ((iterator %in% individuals[, 8]) == FALSE){
      iterator <- iterator + 1
    }
    if ((iterator %in% individuals[, 8]) == TRUE){
      count = length(individuals[, 8][individuals[,8 ] == iterator])
      colony_population <- append(colony_population, count)
      colony_id <- append(colony_id, iterator)
      colony_live <- append(colony_live, 1)
      iterator <- iterator + 1
    }
    else{
      print("there is some other issue you need to figure out")
    }
  }
  colony_numbers <- data.frame(colony_id, colony_population, colony_live)
  #rownames(colony_numbers) <- colony_names
  colnames(colony_numbers) <- c("colony_id", "colony_population", "alive")
  return(colony_numbers)
}


##########################################
## for loop for model run
##########################################

#profvis({
for (i in 1:num_time_steps) {
  colony <- numColonies(individuals)
  num_colonies <- nrow(colony)
  
  ## colony death/dilution
  random <- runif(num_colonies, 0, 100)
  colony[seq_len(num_colonies),][random <= 100*(0.5/24), 3] <- 0 #CMG CHANGED
  
  ## remove individuals from colonies that washed out & update S
  washed_out_colonies <- which(colony[, 3] == 0)
  for (i in 1:num_individuals){
    if ((individuals[i, 8] %in% washed_out_colonies) == TRUE){
      individuals[i, 3] = 0
    }
  }
  
  ## filter by live individuals 
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  individuals_after_death <- length(live_individuals)
  intracellular_after_death <- sum(individuals[live_individuals, 5])
  #individuals_alive <- individuals[!is.na(individuals[,3]) & individuals[,3] == 1,]
  
  ## determine uptake for live individuals
  uptake_per_individual <- vmax*(state$S /(km + state$S))
  foo <- (state$S/(km + state$S))
  ## multiply uptake per individual by each individual's cell size
  resource_uptake_per_individualt <- uptake_per_individual*individuals[live_individuals, 2]
  #resource_uptake_all_ind <- resource_uptake_per_agent*Sr  ### DONT NEED THIS FOR INDIVIDUALS
  ## sum uptake per individual for extracellular mass balance
  uptake_total <- sum(resource_uptake_per_individual)
  
  ## incorporate uptake into individuals
  if (state$S < uptake_total){
    uptake_S <- (state$S/(length(live_individuals)*Sr))
    individuals[live_individuals, 5] <- individuals[live_individuals, 5] + uptake_S
    uptake_total = uptake_S
    state$S = 0 
  }
  else{
    individuals[live_individuals, 5] <- individuals[live_individuals, 5] + resource_uptake_per_agent
    state$S = state$S - uptake_total
  }
  ## update S by dilution
  state$S <- state$S - (state$flow_rate*state$S)
  
  intracellular_after_uptake_before_division <- sum(individuals[live_individuals, 5],na.rm=TRUE)
  
  ## calculate mu (growth rate), convert intracellular resource to quota (nmol/nmol C)
  mu <- individuals[live_individuals, 6]*(1-(individuals[live_individuals, 7]/
                                     (individuals[live_individuals, 5]/individuals[live_individuals, 2])))
  ## update individual cell size of all individuals
  individuals[live_individuals, 2] <- individuals[live_individuals, 2] + (mu*individuals[live_individuals, 2])
  ## if size calculated by division function is below minimum cell size (0.75), 
  ## replace with minimum cell size
  if (any(!is.na(individuals[, 2]) & individuals[, 2] < 0.00022)) { 
    warning("size too small")
  }
  
  ## division by cell size
  
  last_ID <- tail(live_individuals, n = 1)
  for (j in live_individuals){
    if (individuals[j, 2] > (2*m0)){
      individuals[j, 2] = individuals[j, 2]/2 # size update
      
      individuals[j, 5] = individuals[j, 5]/2 # intracellular resource update
      new <- c(last_ID + 1, individuals[j, 2], 1, state$S, individuals[j, 5], individuals[j, 6],
               individuals[j, 7], individuals[j, 8])
      num_individuals <- num_individuals + 1
      total_individuals <- (which(!is.na(individuals[,1])))
      individuals[length(total_individuals) + 1,] <- new
      last_ID <- last_ID + 1
    }
  }
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  individuals_after_division <- length(individuals[live_individuals])
  
  ##update S with new inflow nutrients
  state$S <- state$S + (state$flow_rate*state$nutrient_inflow) 
  
  
  ## update S in individuals matrix
  individuals[, 4] <- state$S
  #individuals_alive <- individuals[!is.na(individuals[,3]) & individuals[,3] == 1,]
  out_per_time <- data.frame("ID" = individuals[live_individuals, 1], 
                             "vmax" = vmax, 
                             "km" = km,
                             "S" = state$S, 
                             "mumax" = individuals[live_individuals, 6],
                             "qnaught" = individuals[live_individuals, 7],
                             "time" = rep(i, nrow(individuals[live_individuals,])), 
                             "intracellular_resource" = individuals[live_individuals, 5],
                             "mu" = mu,
                             "cell_size" = individuals[live_individuals, 2],
                             "alive" = individuals[live_individuals, 3])

  ## filter by live individuals again after division
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  intracellular_resource <- individuals[live_individuals, 5]
  cell_size <- individuals[live_individuals, 2]
  individuals <- matrix(NA, nrow = matrix_size, ncol = 7)
  individuals[1:num_individuals, 1] <- 1:num_individuals  #ID
  individuals[1:num_individuals, 2] <- cell_size[1]  ## cell size
  individuals[1:num_individuals, 3] <- 1 ## alive
  individuals[1:num_individuals, 4] <- state$S ## S
  individuals[1:num_individuals, 5] <- intracellular_resource[1] ## intracellular resource (nmol/cell)
  individuals[1:matrix_size, 6] <- 0.04166667 ## mumax (1/hour)
  individuals[1:matrix_size, 7] <- 0.002333 ## qnaught (nmol/nmol C)
  
  
  ## creating timestep data frames for diagnostics and plotting
  resource_per_time <- data.frame("time" = i, "concentration" = state$S)
  intracellular_per_timestep <- data.frame("time" = i, "intracellular_after_uptake_before_division" = intracellular_after_uptake_before_division)
  intermediate_intracellular <- data.frame("time" = i, "intracellular_before_death" = intracellular_before_death, 
                                           "intracellular_after_death" = intracellular_after_death)
  uptake_per_timestep <- data.frame("time" = i, "uptake_size" = uptake_total)
  uptake_per_individual_timestep <- data.frame("time" = i, "uptake_per_agent" = uptake_per_agent[1])
  Sr_per_timestep <- data.frame("time" = i, "Sr" = Sr, "num_individuals" = num_individuals,
                                "individuals_before_death" = individuals_before_death, "individuals_after_death" = individuals_after_death,
                                "individuals_after_division" = individuals_after_division)
  foo_timestep <- data.frame("time" = i, "S" = state$S, "foo" = foo)
  
  uptake_time <- rbind(uptake_time, uptake_per_timestep)
  resource_time <- rbind(resource_time, resource_per_time)
  intracellular_time <- rbind(intracellular_time, intracellular_per_timestep)
  uptake_per_individual <- rbind(uptake_per_individual, uptake_per_individual_timestep)
  Sr_time <- rbind(Sr_time, Sr_per_timestep)
  intermediate_intracellular_concentration <- rbind(intermediate_intracellular_concentration, intermediate_intracellular)
  foo_time <- rbind(foo_time, foo_timestep)
  
  #row-append to output dataframe that stores outputs for all individuals and all time steps
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

ggplot(Sr_time, aes(x = time, y = (num_individuals))) + geom_point() ##USEFUL
ggplot(Sr_time, aes(x = time, y = (num_individuals*Sr))) + geom_point() ##USEFUL
ggplot(resource_time, aes(x = time, y = concentration)) + geom_point()
ggplot(intracellular_time, aes(x = time, y = intracellular_per_time)) + geom_point(aes(color = "intracellular concentration")) +
  geom_point(aes(y = concentration, color = "S"))
ggplot(intracellular_time, aes(x = time, y = intracellular_after_uptake_before_division)) + geom_point()##USEFUL
ggplot(intracellular_time, aes(x = time, y = mass_balance)) + geom_point()##USEFUL

ggplot(uptake_time, aes(x = time, y = uptake_size)) + geom_point()
ggplot(Sr_time, aes(x = individuals_after_death, y = individuals_after_division)) + geom_point()

ggplot(foo_time, aes(x = S, y = foo)) + geom_point()

plot(output_df$mu ~ output_df$q) #DOES NOT WORK??