

library(GLDEX)
library(ggplot2)
library(dplyr)


##########################################
## Initial individuals and time step parameters
##########################################

num_individuals <- 2000
max_colony_size <- 8000
#Sr = 6666 ## number of individuals a individuals represents
num_time_steps <- 300
vmax <- 0.097083333 #(nmol/nmol C*hour)
km <- 510 #(nmol/L)
m0 <- 0.00022 #(nmole C/cell)
individual_matrix_size <- 1000000

colony_matrix_size <- 1000

## Set up data frame for state variables
state <- data.frame("S" = 10)
state$flow_rate <- 0.02083333
state$nutrient_inflow <- 10

##########################################
## Set up colony for individuals 
##########################################

individuals <- matrix(NA, nrow = individual_matrix_size, ncol = 11)
individuals[1:num_individuals, 1] <- 1:num_individuals  #ID
individuals[1:num_individuals, 2] <- 0.00033  ## cell size
individuals[1:num_individuals, 3] <- 1 ## alive
individuals[1:num_individuals, 4] <- state$S ## S
individuals[1:num_individuals, 5] <- 0.00000079 ## internal resource concentration (nmol/cell)
individuals[1:num_individuals, 6] <- 0.04166667 ## mumax (1/hour)
individuals[1:num_individuals, 7] <- 0.002333 ## qnaught (nmol/nmol C)
individuals[1:num_individuals, 8] <- sample(1:10, num_individuals, replace = TRUE, 
                                                   prob = c(1.5, 0.75, rep(1,8))) ## colony id
individuals[1:num_individuals, 9] <- 0 ## uptake (per individual, 0 to start)
individuals[1:num_individuals, 10] <- 0.097083333 ## Vmax (nmol/nmol C*hour) 
individuals[1:num_individuals, 11] <- 510 ## Km (nmol/L)
  
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
Sr_time <- data.frame()
output_df <- data.frame()

##########################################
## functions created for model
##########################################

numColonies <- function(individuals){
  num_colonies <- length(unique(individuals[1:num_individuals, 8]))
  colony_id <- c()
  colony_population <- c()
  colony_sr <- c()
  for (number in 1:num_colonies){
    count = length(which(!is.na(individuals[,8]) & individuals[,8] == number))
    colony_population <- append(colony_population, count)
    colony_id <- append(colony_id, number)
    colony_sr <- append(colony_sr, 20)
    
  }
  colony_numbers <- data.frame(colony_id, colony_population, colony_sr)
  #rownames(colony_numbers) <- colony_names
  colnames(colony_numbers) <- c("colony_superindividual_id", "colony_population", "Sr")
  return(colony_numbers)
}

##########################################
## for loop for model run
##########################################

#profvis({
for (step in 1:num_time_steps) {
  ## switched this to if statement so Sr isn't reset every timestep
  if (step == 1){
    colony <- numColonies(individuals)
    num_colonies <- nrow(colony)}

  ### 1. DILTUION (of colonies & S) ########################
  for (dilute in 1:(nrow(colony))){
    random <- runif(colony[dilute, 3], 0, 100)
    num_to_remove <- length(which(random <= 100*(0.5/24)))
    colony[dilute, 3] <- (colony[dilute, 3] - num_to_remove)
  }
  # for (remove in 1:nrow(colony)){
  #   if (colony[remove, 3] == 0 & !is.na(colony[remove, 3])){
  #     individuals <- individuals[individuals[,8] != colony[remove, 1], ]
  #     #colony <- colony[colony[, 1] != colony[remove, 1], ]
  #     #num_colonies <- nrow(colony)
  #     #rownames(colony) <- 1:nrow(colony)
  #   }
  #   colony <- colony[colony[, 3] != 0, ]
  #   num_colonies <- nrow(colony)
  # }
  
  remove_list <- colony[colony[,3]==0,1]
  individuals <- subset(individuals,!(individuals[, 8] %in% remove_list))
  colony <- colony[colony[, 3] != 0, ]
  num_colonies <- nrow(colony)
  
  ## get TOTAL individuals (individuals per colony times Sr)
  colony$total_individuals_Sr <- colony$colony_population*colony$Sr
  
  ## update S by dilution
  state$S <- state$S - (state$flow_rate*state$S)
  ######################################################
  
  ### 2. UPTAKE ########################################
  ## filter by live individuals 
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  
  ## NOTE: line 135 allows for variation in Vmax and Km PER INDIVIDUAL
  uptake_per_individual <- vmax*(state$S /(km + state$S)) ## kept this to check that next line works
  individuals[live_individuals, 9] <- individuals[live_individuals, 10]*(state$S /(individuals[live_individuals, 11] + state$S))
  
  ## multiply uptake per individual by each individual's cell size
  resource_uptake_per_individual <- uptake_per_individual*individuals[live_individuals, 2] ## keep this for incorporating S into individuals!!
  individuals[live_individuals, 9] <- resource_uptake_per_individual
  
  ## aggregate uptake per individual to colony level, multiply Sr by uptake per colony to get uptake
  ## per superindividual
  uptake_per_colony <- aggregate(individuals[live_individuals, 9]~individuals[live_individuals, 8], individuals, sum)
  colony$uptake_per_colony <- uptake_per_colony[, 2]
  #colnames(colony) <- c("colony_superindividual_id", "colony_population", "Sr", "uptake_per_superindividal")
  colony$uptake_per_superindividual <- (colony$uptake_per_colony*colony$Sr)
  
  ## sum uptake per individual for extracellular mass balance
  uptake_total <- sum(colony$uptake_per_superindividual)
  
  ## incorporate uptake into individuals
  if (state$S < uptake_total){
    uptake_S <- (state$S/(sum(colony$total_individuals_Sr)))
    individuals[live_individuals, 5] <- individuals[live_individuals, 5] + uptake_S
    #uptake_total = (uptake_S*(sum(colony$total_individuals_Sr))) ## change uptake total to keep track of actual total when desired total exceeds S
    state$S = 0 
  } else {
    individuals[live_individuals, 5] <- individuals[live_individuals, 5] + resource_uptake_per_individual
    state$S = state$S - uptake_total
  }
  ######################################################
  
  ### 3. GROWTH ########################################
  ## calculate mu (growth rate), convert intracellular resource to quota (nmol/nmol C)
  mu <- individuals[live_individuals, 6]*(1-(individuals[live_individuals, 7]/
                                     (individuals[live_individuals, 5]/individuals[live_individuals, 2])))
  ## update individual cell size of all individuals
  individuals[live_individuals, 2] <- individuals[live_individuals, 2] + (mu*individuals[live_individuals, 2])
  ## if size calculated by division is below minimum cell size (m0), replace with m0)
  if (any(!is.na(individuals[, 2]) & individuals[, 2] < 0.00022)) { 
    warning("size too small")
  }
  ######################################################
  
  ### 4. DIVISION ######################################
  last_ID <- tail(live_individuals, n = 1)
  total_individuals <- (which(!is.na(individuals[,1])))
  for (j in live_individuals){
    if (individuals[j, 2] > (2*m0)){
      individuals[j, 2] = individuals[j, 2]/2 # size update
      
      individuals[j, 5] = individuals[j, 5]/2 # intracellular resource update
      new <- c(last_ID + 1, individuals[j, 2], 1, state$S, individuals[j, 5], individuals[j, 6],
               individuals[j, 7], individuals[j, 8], 0, individuals[j, 10], individuals[j, 11])
      num_individuals <- num_individuals + 1
      individuals[length(total_individuals) + 1,] <- new
      total_individuals <- append(total_individuals, length(total_individuals + 1))
      last_ID <- last_ID + 1
    }
  }
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  ### 5. SLOUGHING ################################
  
  for (i in 1:num_colonies){
    if (colony[i, 2] > max_colony_size){
      row.names(colony) <- seq(1, nrow(colony), 1)
      random_sloughing <- runif(colony[i, 2], 0, 100)
      position_to_slough <- which(random_sloughing < 60)
      colony_to_slough <- which(!is.na(individuals[,8]) & individuals[,8] == colony[i, 1])
      sloughed_individual_positions <- colony_to_slough[position_to_slough]
      message("colony sloughed: ", colony[i, 1]," number of sloughed individuals: ", length(sloughed_individual_positions))
      new_colony_id <- max(colony[, 1]) + 1
      individuals[sloughed_individual_positions, 8] = new_colony_id
      if (is.na(sloughed_individual_positions[1])){
        ugh <- sloughed_individual_positions
        print(head(random_sloughing))
        print(head(position_to_slough))
        print(head(colony_to_slough))
        break
      }
      print(head(sloughed_individual_positions))
      new_pop <- length(which(!is.na(individuals[,8]) & individuals[,8] == new_colony_id))
      print(new_pop)
      colony[i, 2] <- colony[i, 2] - new_pop
      new_colony_row <- c(new_colony_id, new_pop, colony[i, 3], (new_pop*colony[i, 3]), 0, 0)
      message("sloughed colony: ", colony[i, 1], " new colony: ", new_colony_id, " time: ", step)
      message(new_colony_row[1], " ", new_colony_row[2], " ", new_colony_row[3], " ", new_colony_row[4])
      print(new_colony_row)
      colony <- rbind(colony, new_colony_row)
      print(tail(colony, n =1))
    }
  }
  # old slough code
  # for (i in 1:num_colonies){
  #   if (colony[i, 2] > max_colony_size){
  #     message("sloughed colony: ", colony[i, 1]," time step: ", step)
  #     random_sloughing <- runif(colony[i, 2], 0, 100)
  #     random_index = 1
  #     while (random_index <= length(random_sloughing)){
  #       for (row in 1:num_individuals){
  #         if ((individuals[row, 8] == colony[i, 1]) && (random_sloughing[random_index] < 60)){
  #           individuals[row, 8] = tail(colony[,1], n=1) + 1
  #           random_index <- random_index + 1
  #           if (random_index >= length(random_sloughing)){
  #             break
  #           }
  #         } else {
  #           random_index <- random_index + 1
  #           if (random_index >= length(random_sloughing)){
  #             break
  #           }
  #         }
  #       }
  #       message("sloughed colony id: ", colony[i,1], " sloughed colony pop after: ", length(which(!is.na(individuals[,8]) & individuals[,8] == colony[i,1])))
  #       message("new colony id:", tail(colony[,1], n=1) + 1, " new colony pop: ", length(which(!is.na(individuals[,8]) & individuals[,8] == tail(colony[,1], n=1) + 1)))
  #     }
  #     new_pop <- length(which(!is.na(individuals[,8]) & individuals[,8] == tail(colony[,1], n=1) + 1))
  #     new_colony_row <- c(tail(colony[,1], n=1) + 1, new_pop, colony[i, 3], (new_pop*colony[i, 3]), 0, 0)
  #     colony <- rbind(colony, new_colony_row)
  #   }
  # }
  
  live_individuals <- which(!is.na(individuals[,3]) & individuals[,3] == 1)
  #individuals_after_division <- length(individuals[live_individuals])
  
  ## update number of colonies
  num_colonies <- nrow(colony)
  
  ## create vector of updated colony populations
  new_colonies <- tabulate(individuals[live_individuals, 8])
  new_colonies <- fun.zero.omit(new_colonies)
  
  #print(tail(colony,n=1))
  
  ## update colony populations
  for (k in 1:num_colonies){
    colony[k, 2] <- new_colonies[k,]
  }
  
  ## another way to subset
  ##colony[colony$… == …]
  ##colony$colony_pop[colony$id==12]<- ...
  
  #print(tail(colony,n=1))
  
  ######################################################
  
  ### 6. MERGING/SEPARATING ############################
  ######################################################
  
  ### 7. UPDATE S W/ NEW RESOURCE ######################
  state$S <- state$S + (state$flow_rate*state$nutrient_inflow) 
  
  ## update S in individuals matrix
  individuals[, 4] <- state$S
  ######################################################
  
  ### 8. OUTPUT DATAFRAME FOR ALL TIMESTEPS ############
  out_per_time <- data.frame("ID" = individuals[live_individuals, 1],
                             "time" = rep(step, nrow(individuals[live_individuals,])),
                             "colony_ID" = individuals[live_individuals, 8],
                             "S" = state$S,
                             "intracellular_resource" = individuals[live_individuals, 5],
                             "cell_size" = individuals[live_individuals, 2],
                             "vmax" = individuals[live_individuals, 10], 
                             "km" = individuals[live_individuals, 11],
                             "mumax" = individuals[live_individuals, 6],
                             "qnaught" = individuals[live_individuals, 7],
                             "mu" = mu)

  ## recreate individuals matrix for new timestep ### GET RID OF THIS????
  # intracellular_resource <- individuals[live_individuals, 5]
  # cell_size <- individuals[live_individuals, 2]
  # individuals <- matrix(NA, nrow = matrix_size, ncol = 7)
  # individuals[1:num_individuals, 1] <- 1:num_individuals  #ID
  # individuals[1:num_individuals, 2] <- cell_size[1]  ## cell size
  # individuals[1:num_individuals, 3] <- 1 ## alive
  # individuals[1:num_individuals, 4] <- state$S ## S
  # individuals[1:num_individuals, 5] <- intracellular_resource[1] ## intracellular resource (nmol/cell)
  # individuals[1:matrix_size, 6] <- 0.04166667 ## mumax (1/hour)
  # individuals[1:matrix_size, 7] <- 0.002333 ## qnaught (nmol/nmol C)
  
  
  ## creating timestep data frames for diagnostics and plotting
  resource_per_time <- data.frame("time" = step, "concentration" = state$S)

  uptake_per_timestep <- data.frame("time" = step, "uptake_size" = uptake_total)
  
  uptake_time <- rbind(uptake_time, uptake_per_timestep)
  resource_time <- rbind(resource_time, resource_per_time)
  
  #row-append to output dataframe that stores outputs for all individuals and all time steps
  output[[i]] <- out_per_time
  output_df <- rbind(output_df, out_per_time)
  ######################################################
}
#})

##########################################
## Plotting
##########################################

output_df <- data.frame()

for (i in 1:num_time_steps){
  output_interem <- as.data.frame(output[[i]])
  output_df <- rbind(output_df, output[[i]])
}

intracellular_time <- cbind(intracellular_time, resource_time[,2, drop=FALSE])
intracellular_time$mass_balance <- (intracellular_time$intracellular_after_uptake_before_division+intracellular_time$concentration)


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

plot(output_df$cell_size ~ output_df$time)
plot(output_df$intracellular_resource ~ output_df$time)
plot(output_df$S ~ output_df$time)

plot(output_df$mu ~ output_df$q) #DOES NOT WORK??