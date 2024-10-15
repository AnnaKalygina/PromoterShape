library(ggplot2)
library(reshape2)
library(gridExtra)

MGW <- read.table("~/Downloads/athaliana_MGW_200.txt", header = FALSE)

# Calculate the average of each column
average_data <- colMeans(data)

positions <- seq(-200, 200, length.out = ncol(data))

average_df <- data.frame(Position = positions, Mean_Value = average_data)

library(ggplot2)
ggplot(average_df, aes(x = Position, y = Mean_Value)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(x = "Position around TSS (bp)", y = "Average MGW", title = "Minor Groove Width") +
  scale_x_continuous(breaks = seq(-200, 200, 50)) 




properties <- c("MGW", "Buckle", 
                "Opening", "Tilt")


input_dir <- "~/Downloads/athaliana_"
plot_list <- list()

for (prop in properties) {
  
  # Construct file path for each property
  file_path <- paste0(input_dir, prop, "_200.txt")
  
  # Read the data
  data <- read.table(file_path, header = FALSE)
  
  # Calculate the average of each column
  average_data <- colMeans(data)
  
  # Create a new Position vector ranging from -200 to 200
  positions <- seq(-200, 200, length.out = ncol(data))
  
  # Create a data frame for plotting
  average_df <- data.frame(Position = positions, Mean_Value = average_data)
  
  # Plot the averaged data
  p <- ggplot(average_df, aes(x = Position, y = Mean_Value)) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(x = "Position", y = "Average Feature Value", 
         title = paste("Average DNA Shape Prediction for", prop)) +
    scale_x_continuous(breaks = seq(-200, 200, 50))  # Adjust x-axis ticks
  
  # Store the plot in the list
  plot_list[[prop]] <- p
}

# Arrange the plots in a 1-column and 4-row layout using gridExtra
grid.arrange(grobs = plot_list, ncol = 1)
