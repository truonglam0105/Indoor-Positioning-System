#_______________________________________________________________________________#
#_______________________________________________________________________________#
## LOADING, FORMATTING, AND CLEANING THE OFFLINE DATA:
#_______________________________________________________________________________#
#_______________________________________________________________________________#
# Install once, only if packages are not already installed
#install.packages("dplyr")

# Load libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(magrittr)

#-------------------------------------------------------------------------------#
# sets the current working directory
setwd("C:/Users/risin/OneDrive/Desktop/School/Classes/STAT410/IPS_data")
#-------------------------------------------------------------------------------#
# Loads the offline.final.trace.txt into an object
offline_data <- readLines("offline.final.trace.txt")
#-------------------------------------------------------------------------------#
# offline_data is a character vector of length 151392
str(offline_data)
#-------------------------------------------------------------------------------#
# There are 5312 comments in the data
sum(substr(offline_data,1,1) =="#")
#-------------------------------------------------------------------------------#
# This function splits the text line to match specific columns names its respective row of data
processLine <- function(x){
  # Regex in the strsplit function to split data on ;=, characeters (from looking at the data)
  tokens = strsplit(x,"[;=,]")[[1]]
  # Handles the case where there are no recorded signals (not useful for this project) and removes the data                           
  if (length(tokens) == 10){
    return(NULL)
  }
  # For each signal recording, tokens 1,3,5, and 9 are columns names
  tmp = matrix(tokens[-(1:10)],ncol=4, byrow=T)
  cbind(matrix(tokens[c(2,4,6:8,10)],nrow=nrow(tmp),ncol=6,byrow=T),tmp)
}
#-------------------------------------------------------------------------------#
# Runs the porcessLine function over the entire data set to build a data frame
lines <- offline_data[substr(offline_data,1,1) != "#"] # Removes comments from data
tmp <- lapply(lines,processLine) 
train <- as.data.frame(do.call("rbind",tmp), stringAsFactors=F)
#-------------------------------------------------------------------------------#
# Check dimensions of training set
dim(train)
head(train)
str(train)
#-------------------------------------------------------------------------------#
# Assigns names to each of the columns
names(train) <- c("time","scanMac","posX","posY","posZ","orientation","mac","signal","channel","type")
#-------------------------------------------------------------------------------#
# Checks for any NA's in the data set
anyNA(train)
#-------------------------------------------------------------------------------#
# Changes each variable to the correct type
train_type <- train %>%
  mutate(across(c(time, posX, posY, posZ, orientation, signal, channel, type), ~ suppressWarnings(as.numeric(.))))
#-------------------------------------------------------------------------------#
# Removes variables that provide redundant or no information
train_no_scanMac <- train_type %>% select(-scanMac,-posZ)
#-------------------------------------------------------------------------------#
# Creates a function to convert the observations in the time variable from milliseconds to into POSIXct time format 
posixct_func<- function(x){  
  seconds<- x %/% 1000
  date_time<-as.POSIXct(seconds, origin = "1970-01-01", tz = "GMT")
}
#-------------------------------------------------------------------------------#
# Creates a new data frame with with time converted into posixct format
train_posixct <- train_no_scanMac %>% mutate(time = posixct_func(time))
#-------------------------------------------------------------------------------#
# Converts the accessPointLocations.txt file information into a table  
mac_locs <- readr::read_table("accessPointLocations.txt")
#-------------------------------------------------------------------------------#
# Stores the information from the Macs column in the mac_locs table
MAC_addresses<-c(mac_locs$Macs)
#-------------------------------------------------------------------------------#
# Filters the data for instances of relevant access points only
train_mac <- train_posixct %>% filter(mac %in% MAC_addresses)
#-------------------------------------------------------------------------------#
# renames the Macs column to match the data set
colnames(mac_locs)[colnames(mac_locs) == "Macs"] <- "mac"
#-------------------------------------------------------------------------------#
# use inner join to find the matched data effectively adding x and y columns 
train_access_points <- train_mac %>% 
  inner_join(mac_locs, by = c("mac")) %>% 
  select(
    1:5,          # Keeps the first 5 columns in place
    macX = x,     # Moves/renames the mac x coordinates to column 6 
    macY = y,     # Moves/renames the mac y coordinates to column 7
    everything()  # Moves everything else in order after that
  )
#-------------------------------------------------------------------------------#
# This function takes and an angle and locates its closest proximity to angles in the bin vector
nearest_angle <- function(angle, bins) {
  
  selected_val <- bins[which.min(abs(angle - bins))]
  if(selected_val != 360){
    return(selected_val)
  } else{
    return(0)
  }
}
bins <- c(0, 45, 90, 135, 180, 225, 270, 315, 360)
#-------------------------------------------------------------------------------#
# Applies the nearest_angle function to the orientation column and creates a group version of the orientation variable conformed to angles 0, 45, 90, 135, 180, 225, 270, 315.
train_orientation <- train_access_points %>%
  transmute(
    across(1:3),                  # Retain the first three columns
    adj_orient = sapply(orientation, nearest_angle, bins = bins),  # Create the new column
    across(-orientation)          # Retain all other columns except 'orientation'
  )
#-------------------------------------------------------------------------------#
# Find the Euclidean distance between location of the device and access point
dist <- sqrt((train_orientation$posX - train_orientation$macX)^2 + (train_orientation$posY - train_orientation$macY)^2)
dist <- round(dist, digits = 2)
#-------------------------------------------------------------------------------#
# Add 'dist' column after the 6th column
train_dist <- train_orientation %>%
  add_column(dist = dist, .after = 7)
#-------------------------------------------------------------------------------#
# Consolidates useful information by groups and performs summary statistics   
train_summary <- train_dist %>%
  group_by(posX, posY, adj_orient, mac, macX, macY, dist) %>%
  summarise(mean_signal = mean(signal, na.rm = TRUE),
            median_signal = median(signal, na.rm = TRUE),
            sd_signal = sd(signal, na.rm = TRUE),
            min_signal = min(signal, na.rm = TRUE),
            max_signal = max(signal, na.rm = TRUE),
            IQR_signal = IQR(signal, na.rm = TRUE),
            Q1 = quantile(signal, 0.25, na.rm = TRUE),        
            Q3 = quantile(signal, 0.75, na.rm = TRUE),
            IQR_min = Q1 - 1.5 * IQR_signal,          
            IQR_max = Q3 + 1.5 * IQR_signal,
            outliers = sum(signal < IQR_min |      # Counts outliers
                             signal > IQR_max),
            count = n(),
            .groups = "drop"  # Removes grouping from the resulting data frame
  )%>% # Removes information no longer needed after getting the ouliers
  select(-min_signal,-max_signal,-IQR_signal,-Q1,-Q3,-IQR_min,-IQR_max)
#head(train_summary)
#-------------------------------------------------------------------------------#
# Calculate the angle in radians from device to router
angle <- atan2(train_summary$macY - train_summary$posY, train_summary$macX - train_summary$posX)
# If desired, convert the angle to degrees
angle <- angle * (180 / pi)
# Adjust to ensure all angles are in the range 0-359 degrees
angle <- (angle + 360) %% 360
# Round to the 6th decimal place
angle <- round(angle, digits = 2)
train_angle <- data.frame(append(train_summary, list(angle = angle), after = 3))
#head(train_angle,48)
#-------------------------------------------------------------------------------#
train_angle_diff <- train_angle %>%
  group_by(posX, posY, mac, adj_orient) %>%
  # Step 2: Calculate the angular difference, considering the circular nature
  mutate(
    angle_diff = abs(adj_orient - angle),
    # Correct for values greater than 180 degrees by subtracting from 360
    angle_diff = ifelse(angle_diff > 180, 360 - angle_diff, angle_diff)) %>% 
  select(
    1:4,           # select the first 4 columns to be the first
    angle_diff,    # moves angle_diff to the 5th column spot      
    everything()   # Then adds the rest of the columns in order
  )
#print(train_angle_diff, n = 30)

# Checks to make sure there are no angle differences below zero or above 180, or NA values 
unique_angle_diff<-unique(train_angle_diff$angle_diff)
# Logical condition for values below 0, above 180, or NA
invalid_cases <- unique_angle_diff < 0 | unique_angle_diff > 180 | is.na(unique_angle_diff)
# Extract the problematic values
problem_values <- unique_angle_diff[invalid_cases]
#print(problem_values)
#-------------------------------------------------------------------------------#
# Finds the unique combinations of locations and mac addressees 
# where the device is pointing in the direction closest to the router/access point
train_min_adiff <- train_angle_diff %>% 
  group_by(posX,posY,mac) %>% 
  summarise(
  angle_diff = min(angle_diff),  # Find the minimum value for the angle difference
  .groups = "drop"  # grouping is not necessary at this point
)
# Matches data from test_min_adiff and test_angle_diff so we have the data
# we want all the columns where the device was pointing closes to the router in each recording
train_improved_sig_accuracy<- train_angle_diff %>% semi_join(train_min_adiff)
#print(train_improved_sig_accuracy, n = 18)
#-------------------------------------------------------------------------------#
# Summarizes the 110 measurements for each location/orientation group to be able use the median signal
train_closestAPs_any_orient <- train_summary %>%
  group_by(posX, posY, dist) 
#-------------------------------------------------------------------------------#
# Creates a new data frames showing the 3 closest access points associated with each location (where the device is facing the router)
train_closestAPs_facing_router <- train_improved_sig_accuracy %>%
  group_by(posX, posY) %>%                      # Regroup by position
  arrange(dist, .by_group = TRUE) %>%           # Sort by distance
  slice(1:3)     
#print(train_closestAPs_facing_router, n = 18)
#-------------------------------------------------------------------------------#
# Creates a new data frames showing the 3 closest access points associated with each location (where the signal strength is the best)
train_closestAPs_max_signal <- train_summary %>%
  group_by(posX, posY, mac) %>%                 # Group by position and router
  filter(mean_signal == max(mean_signal)) %>%   # Keep the recording with the strongest signal per router
  ungroup() %>%                                 # Remove grouping to avoid interference
  group_by(posX, posY) %>%                      # Regroup by position
  arrange(dist, .by_group = TRUE) %>%           # Sort by distance
  distinct(mac, .keep_all = TRUE) %>%           # Keep only one row per unique router
  slice(1:3)      
#-------------------------------------------------------------------------------#
# Rewrites the new information over the original testing data for ease of code
train_final1 <- train_dist
train_final2 <- train_closestAPs_any_orient
train_final3 <- train_closestAPs_facing_router
train_final4 <- train_closestAPs_max_signal
#head(train_final1,n = 5, width = Inf)
#head(train_final2,n = 5,width = Inf)
#head(train_final3,n = 5,width = Inf)
#head(train_final4,n = 5,width = Inf)
#-------------------------------------------------------------------------------#
# Plot the relationship between signal strength and distance for the entire data set
plot(train_final1$signal, train_final1$dist, main = "Offline Distance vs Signal \n (all data)", xlab = "Offline_Signal", ylab = "Offline_Distance", pch = 16)
#-------------------------------------------------------------------------------#
# Plot the relationship between signal strength and distance when noise is removed (for all orientations)
plot(train_final2$median_signal, train_final2$dist, main = "Offline Distance vs Signal \n (all six routers with all orientions)", xlab = "Offline_Signal", ylab = "Offline_Distance", pch = 16)
#-------------------------------------------------------------------------------#
# Plot the relationship between signal strength and distance when noise is removed (and the when device is pointed at the router)
plot(train_final3$median_signal, train_final3$dist, main = "Offline Distance vs Signal \n (3 closest routers with device pointed towards router)", xlab = "Offline_Signal", ylab = "Offline_Distance", pch = 16)
#-------------------------------------------------------------------------------#
# Plot the relationship between signal strength and distance when noise is removed (and when the signal strength is the best) 
plot(train_final4$median_signal, train_final4$dist, main = "Offline Distance vs Signal \n (3 closest routers with best signal strength)", xlab = "Offline_Signal", ylab = "Offline_Distance", pch = 16)
#_______________________________________________________________________________#
#_______________________________________________________________________________#
## LOADING, FORMATTING, AND CLEANING THE ONLINE DATA:
#_______________________________________________________________________________#
#_______________________________________________________________________________#
# Loads the online.final.trace.txt into an object
online_data <- readLines("online.final.trace.txt")
#-------------------------------------------------------------------------------#
# online_data is a character vector of length 6832
str(online_data)
#-------------------------------------------------------------------------------#
# There are 240 comments in the data
sum(substr(online_data,1,1) =="#")
#-------------------------------------------------------------------------------#
# Runs the porcessLine function over the entire data set to build a data frame
lines <- online_data[substr(online_data,1,1) != "#"] # Removes comments from data
tmp <- lapply(lines,processLine) 
test <- as.data.frame(do.call("rbind",tmp), stringAsFactors=F)
#-------------------------------------------------------------------------------#
# Check dimensions of training set
dim(test)
head(test)
str(test)
#-------------------------------------------------------------------------------#
# Assigns names to each of the columns
names(test) <- c("time","scanMac","posX","posY","posZ","orientation","mac","signal","channel","type")
#-------------------------------------------------------------------------------#
# Checks for any NA's in the data set
anyNA(test)
#-------------------------------------------------------------------------------#
# Changes each variable to the correct type
test_type <- test %>%
  mutate(across(c(time, posX, posY, posZ, orientation, signal, channel, type), ~ suppressWarnings(as.numeric(.))))
#-------------------------------------------------------------------------------#
# Removes variables that provide redundant or no information
test_no_scanMac <- test_type %>% select(-scanMac,-posZ)
#-------------------------------------------------------------------------------#
# Creates a function to convert the observations in the time variable from milliseconds to into POSIXct time format 
posixct_func<- function(x){  
  seconds<- x %/% 1000
  date_time<-as.POSIXct(seconds, origin = "1970-01-01", tz = "GMT")
}
# Creates a new data frame with with time converted into posixct format
test_posixct <- test_no_scanMac %>% mutate(time = posixct_func(time))
#-------------------------------------------------------------------------------#
# Converts the accessPointLocations.txt file information into a table  
mac_locs <- readr::read_table("accessPointLocations.txt")
#-------------------------------------------------------------------------------#
# Stores the information from the Macs column in the mac_locs table
MAC_addresses<-c(mac_locs$Macs)
#-------------------------------------------------------------------------------#
# Filters the data for instances of relevant access points only
test_mac <- test_posixct %>% filter(mac %in% MAC_addresses)
#-------------------------------------------------------------------------------#
# renames the Macs column to match the data set
colnames(mac_locs)[colnames(mac_locs) == "Macs"] <- "mac"
#-------------------------------------------------------------------------------#
# use inner join to find the matched data effectively adding x and y columns 
test_access_points <- test_mac %>% 
  inner_join(mac_locs, by = c("mac")) %>% 
  select(
    1:5,          # Keeps the first 5 columns in place
    macX = x,     # Moves/renames the mac x coordinates to column 6 
    macY = y,     # Moves/renames the mac y coordinates to column 7
    everything()  # Moves everything else in order after that
  )
#-------------------------------------------------------------------------------#
# This function takes and an angle and locates its closest proximity to angles in the bin vector
nearest_angle <- function(angle, bins) {
  
  selected_val <- bins[which.min(abs(angle - bins))]
  if(selected_val != 360){
    return(selected_val)
  } else{
    return(0)
  }
}
bins <- c(0, 45, 90, 135, 180, 225, 270, 315, 360)
#-------------------------------------------------------------------------------#
# Applies the nearest_angle function to the orientation column and creates a group version of the orientation variable conformed to angles 0, 45, 90, 135, 180, 225, 270, 315.
test_orientation <- test_access_points %>%
  transmute(
    across(1:3),                  # Retain the first three columns
    adj_orient = sapply(orientation, nearest_angle, bins = bins),  # Create the new column
    across(-orientation)          # Retain all other columns except 'orientation'
  )
#-------------------------------------------------------------------------------#
# Find the Euclidean distance between location of the device and access point
dist <- sqrt((test_orientation$posX - test_orientation$macX)^2 + (test_orientation$posY - test_orientation$macY)^2)
dist <- round(dist, digits = 2)
#-------------------------------------------------------------------------------#
# Add 'dist' column after the 6th column
test_dist <- test_orientation %>%
  add_column(dist = dist, .after = 7)
#-------------------------------------------------------------------------------#
# Consolidates useful information by groups and performs summary statistics   
test_summary <- test_dist %>%
  group_by(posX, posY, adj_orient, mac, macX, macY, dist) %>%
  summarise(mean_signal = mean(signal, na.rm = TRUE),
            median_signal = median(signal, na.rm = TRUE),
            sd_signal = sd(signal, na.rm = TRUE),
            min_signal = min(signal, na.rm = TRUE),
            max_signal = max(signal, na.rm = TRUE),
            IQR_signal = IQR(signal, na.rm = TRUE),
            Q1 = quantile(signal, 0.25, na.rm = TRUE),        
            Q3 = quantile(signal, 0.75, na.rm = TRUE),
            IQR_min = Q1 - 1.5 * IQR_signal,          
            IQR_max = Q3 + 1.5 * IQR_signal,
            outliers = sum(signal < IQR_min |      # Counts outliers
                             signal > IQR_max),
            count = n(),
            .groups = "drop"  # Removes grouping from the resulting data frame
  )%>% # Removes information no longer needed after getting the ouliers
  select(-min_signal,-max_signal,-IQR_signal,-Q1,-Q3,-IQR_min,-IQR_max)
#head(test_summary)
#-------------------------------------------------------------------------------#
# For the three closest routers to each point find the highest median signal value. 
test_closestAPs_max_signal <- test_summary %>%
  group_by(posX, posY, mac) %>%                 # Group by position and router
  filter(mean_signal == max(mean_signal)) %>%   # Keep the recording with the strongest medain signal per router
  ungroup() %>%                                 # Remove grouping to avoid interference
  group_by(posX, posY) %>%                      # Regroup by position
  arrange(dist, .by_group = TRUE) %>%           # Sort by distance
  distinct(mac, .keep_all = TRUE) %>%           # Keep only one row per unique router
  slice(1:3)      
#-------------------------------------------------------------------------------#
# Rewrites the new information over the original testing data for ease of code
test_final <- test_closestAPs_max_signal
#head(test_final,n = 50, width = Inf)
#_______________________________________________________________________________#
#_______________________________________________________________________________#
## PREDICTION MODELING
#_______________________________________________________________________________#
#_______________________________________________________________________________#
## WE DON'T NEED TO NORMALIZE THE PREDICTOR VARIABLES FOR KNN REGRESSION 
## (THERE IS ONLY ONE PREDICTOR VARIABLE FOR THE MODEL)
#-------------------------------------------------------------------------------#
### USING TRAIN_FINAL 1,2, OR 3 YOU WILL SEE DIFFERENT RESULTS (IN KNN REGRESSION)
#-------------------------------------------------------------------------------#
# Load libraries
# install.packages("FNN")
# install.packages("caret")

library(FNN)
library(caret)

# Define a function for KNN regression and residuals analysis
knn_residuals_analysis <- function(train_x, test_x, train_y, test_y, k) {
  # Ensure inputs are data frames
  train_x <- as.data.frame(train_x)
  test_x <- as.data.frame(test_x)
  
  # Train the KNN model
  knn_model <- knn.reg(train = train_x, test = test_x, y = train_y, k = k)
  
  # Get predictions
  predictions <- knn_model$pred
  
  # Compute residuals
  residuals <- test_y - predictions
  
  # Calculate RMSE
  rmse <- RMSE(predictions, test_y)
  
  # Print RMSE
  print(paste("RMSE for KNN model:", rmse))
  
  # Display residuals in a data frame
  residuals_df <- data.frame(
    Observation = seq_along(residuals),
    Actual = test_y,
    Predicted = predictions,
    Residuals = residuals
  )
  print("Residuals Data Frame:")
  print(residuals_df)
  
  # Plot residuals (Vertical Scatter Plot)
  plot(
    residuals,                     # Residuals are the y-values
    xlab = "Observation Index",    # Label for the x-axis
    ylab = "Residuals",            # Label for the y-axis
    main = "Residuals Plot",       # Title for the plot
    pch = 19,                      # Use solid circles for points
    col = "blue"                   # Color of the points
  )
  
  # Create Q-Q plot for residuals
  qqnorm(residuals, 
         main = "Q-Q Plot of Residuals", 
         xlab = "Theoretical Quantiles", 
         ylab = "Sample Quantiles", 
         col = "lightblue")
  qqline(residuals, col = "red", lwd = 2)  # Add reference line in red
  
  # Residuals vs Predicted Values Plot
  plot(
    predictions, residuals,
    xlab = "Predicted Values",
    ylab = "Residuals",
    main = "Residuals vs Predicted Values",
    pch = 19,
    col = "red"
  )
  abline(h = 0, col = "blue", lwd = 2)  # Add horizontal line at y = 0
  
  # Return residuals data frame and RMSE as output
  return(knn_model)
}

# Example usage
# Assuming train_final, test_final, and k are defined
### NOTE!!! 
### Use signal for train_final1 and train_final2
### Use median_signal for train_final3
train_x <- train_final3$median_signal
test_x <- test_final$median_signal
train_y <- train_final3$dist
test_y <- test_final$dist
k <- 19

knn_model <- knn_residuals_analysis(train_x, test_x, train_y, test_y, k)
#-------------------------------------------------------------------------------#

## NOTE: MASS AND DYPLR HAVE CONFLICTS IN THEIR PACKAGES (ESP USING SELECTED FUNCTION)
#install.packages("MASS") # for ginv
#install.packages("pracma") # for pinv

#library(MASS) # for ginv
library(pracma) # for pinv

test_final$pred_dist <- knn_model$pred

test_final <- test_final %>%
  select(
    1:8,
    pred_dist,
    everything()
  )

# Function to compute least squares solution (pinv and ginv when A is singular)
compute_least_squares <- function(data) {
  results <- list()  # To store results for each location
  
  # Group the data by posX and posY
  grouped_data <- split(data, list(test_final$posX, test_final$posY))
  
  for (group in grouped_data) {
    if (nrow(group) == 3) {  # Ensure there are exactly 3 rows per location
      # Extract macX, macY, and distances
      macX <- group$macX
      macY <- group$macY
      d <- group$pred_dist
      
      # Extract reference point (first router)
      x1 <- macX[1]
      y1 <- macY[1]
      d1 <- d[1]
      
      # Extract second and third points
      x2 <- macX[2]
      y2 <- macY[2]
      d2 <- d[2]
      
      x3 <- macX[3]
      y3 <- macY[3]
      d3 <- d[3]
      
      # Construct A matrix
      A <- matrix(c(
        2 * (x2 - x1), 2 * (y2 - y1),
        2 * (x3 - x1), 2 * (y3 - y1)
      ), nrow = 2, byrow = TRUE)
      
      # Construct b vector
      b <- c(
        (d1^2 - d2^2) + (x2^2 - x1^2) + (y2^2 - y1^2),
        (d1^2 - d3^2) + (x3^2 - x1^2) + (y3^2 - y1^2)
      )
      
      # Solve Ax = b using the normal equations: x = (A^T A)^(-1) A^T b
      #AtA_inv <- solve(t(A) %*% A)  # Compute (A^T A)^(-1)
      #Atb <- t(A) %*% b  # Compute A^T b
      #x <- AtA_inv %*% Atb  # Solve for x
      
      # Solve Ax = b using pseudoinverse
      x <- pinv(A) %*% b
      # x <- ginv(A) %*% b
      
      # Save the computed coordinates
      results <- append(results, list(list(posX = group$posX[1], posY = group$posY[1], predX = x[1], predY = x[2])))
    
      } 
  }
  # Convert results to a data frame
  results_df <- do.call(rbind, lapply(results, as.data.frame))
  
  # Calculate the Euclidean distance between the predicted and actual points
  results_df$dist <- sqrt((results_df$posX - results_df$predX)^2 + (results_df$posY - results_df$predY)^2)
  results_df$dist <- round(results_df$dist, digits = 2)
  
  return(results_df)
}

# Function to calculate RMSE
compute_rmse <- function(results_df) {
  # Calculate the squared differences
  squared_diff <- (results_df$predX - results_df$posX)^2 + (results_df$predY - results_df$posY)^2
  
  # Compute the mean of squared differences
  mse <- mean(squared_diff)
  
  # Take the square root to get RMSE
  rmse <- sqrt(mse)
  
  return(rmse)
}

# Get the location predictions of the least squares model
results_df <- compute_least_squares(test_final)

# Apply the function on your results_df
rmse <- compute_rmse(results_df)

# Print the RMSE
print(paste("The RMSE of the least squares method is:", round(rmse, 2)))

# 1. Histogram of Prediction Errors
# 3. Prediction Error Heatmap
ggplot(results_df, aes(x = posX, y = posY, color = dist)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Prediction Error Heatmap",
       x = "X Coordinate", y = "Y Coordinate", color = "Error") +
  theme_minimal()


# Step 1: Calculate residuals
results_df$residualX <- results_df$posX - results_df$predX
results_df$residualY <- results_df$posY - results_df$predY

# Q-Q Plot of Residuals (X)
qqnorm(results_df$residualX, main = "Q-Q Plot of Residuals (X)")
qqline(results_df$residualX, col = "red")

# Q-Q Plot of Residuals (Y)
qqnorm(results_df$residualY, main = "Q-Q Plot of Residuals (Y)")
qqline(results_df$residualY, col = "red")

# Histogram of Residuals for X
ggplot(results_df, aes(x = residualX)) +
  geom_histogram(color = "black", fill = "lightblue", bins = 20) +
  labs(title = "Histogram of Residuals (X)", x = "Residual X", y = "Count") +
  theme_minimal()

# Histogram of Residuals for Y
ggplot(results_df, aes(x = residualY)) +
  geom_histogram(color = "black", fill = "lightgreen", bins = 20) +
  labs(title = "Histogram of Residuals (Y)", x = "Residual Y", y = "Count") +
  theme_minimal()

# 2. Predicted vs Actual Positions with Error Arrows
ggplot(results_df) +
  # Actual positions
  geom_point(aes(x = posX, y = posY, color = "Actual Position"), linewidth = 3, alpha = 0.6) + 
  # Predicted positions
  geom_point(aes(x = predX, y = predY, color = "Predicted Position"), linewidth = 3, alpha = 0.6) + 
  # Arrows between actual and predicted
  geom_segment(aes(x = posX, y = posY, xend = predX, yend = predY), 
                arrow = arrow(length = unit(0.2, "cm")), color = "gray") +
  labs(
    title = "Predicted vs Actual Positions with Error Arrows",
    x = "X Coordinate", 
    y = "Y Coordinate",
    color = "Legend"  # Title for the legend
  ) +
  theme_minimal()
  
# 3. Prediction Error Heatmap
ggplot(results_df, aes(x = posX, y = posY, color = dist)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Prediction Error Heatmap",
        x = "X Coordinate", y = "Y Coordinate", color = "Error") +
  theme_minimal()
  

