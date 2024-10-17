library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)

MGW <- read.table("~/Downloads/athaliana_MGW_200.txt", header = FALSE)

# Collapsing the data
MGW_z_score <- scale(MGW)
MGW_z_score_average <- colMeans(MGW_z_score)
str(MGW_z_score)
col_sds <- apply(MGW_z_score, 2, sd)
print(col_sds)
# Calculate the average of each column
average_data <- colMeans(MGW)

par(mfrow=c(1, 2))
plot(MGW_z_scored, col = "red")
plot(average_data, col = "blue")

positions <- seq(-200, 200, length.out = ncol(data))

average_df <- data.frame(Position = positions, Mean_Value = average_data)

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



MGW_data <- read.table("~/Downloads/athaliana_MGW_200.txt", header = FALSE)

col_means <- colMeans(MGW_data)
col_sds <- apply(MGW_data, 2, sd)

MGW_z_score_manual <- sweep(MGW_data, 2, col_means, FUN = "-") 
MGW_z_score_manual <- sweep(MGW_z_score_manual, 2, col_sds, FUN = "/")  
MGW_z_scored = colMeans(MGW_z_score_manual)


MGW_z_score_rows <- t(apply(MGW_data, 1, function(x) (x - mean(x)) / sd(x)))
MGW_z_scored = colMeans(MGW_z_score_rows)
average_data = colMeans(MGW_data)[100:300]

par(mfrow=c(1, 2))
plot(MGW_z_scored, col = "red")
plot(average_data, col = "blue")

simple <- data.frame(
  Position = seq(-100, 100, length.out = length(avg_z_scores)),  
  Avg_Z_Score = average_data
)

avg_z_scores <- colMeans(MGW_z_score_rows)[100:300]
df <- data.frame(
  Position = seq(-100, 100, length.out = length(avg_z_scores)),  
  Avg_Z_Score = avg_z_scores          
)

a <- ggplot(df, aes(x = Position, y = Avg_Z_Score)) +
  geom_point(color = "blue") +  
  labs(x = "Position around TSS (bp)", y = "Average Z-score") +
  ggtitle("Normalised Minor Groove Width") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
  

b <- ggplot(simple, aes(x = Position, y = Avg_Z_Score)) +
  geom_point(color = "red") +  
  labs(x = "Position around TSS (bp)", y = "Average score") +
  ggtitle(" Not normalised Minor Groove Width") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(a, b, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
