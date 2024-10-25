library(ggplot2)
library(reshape2)
library(gridExtra)

properties <- c("Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt")



input_dir <- "~/Desktop/Promoters/dna_shape/hsapiens_"
plot_list_hsapiens <- list()

for (prop in properties) {
  

  file_path <- paste0(input_dir, prop, "_200.txt")
  data <- read.table(file_path, header = FALSE)
  
  average_data <- colMeans(data)
  
  positions <- seq(-200, 200, length.out = ncol(data))
  
  average_df <- data.frame(Position = positions, Mean_Value = average_data)
  
  p <- ggplot(average_df, aes(x = Position, y = Mean_Value)) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(x = "Position", y = "Average Feature Value", 
         title = paste("Average DNA Shape Prediction for", prop, " in H.sapiens")) +
    scale_x_continuous(breaks = seq(-200, 200, 50))  
  plot_list_hsapiens[[prop]] <- p
}


grid.arrange(grobs = plot_list_hsapiens[1:5], ncol = 1)


### Can we now plot controls?

input_dir <- "~/Desktop/Promoters/dna_shape_control/hsapiens_"
plot_list_hsapiens_hmm <- list()

for (prop in properties) {
  
  
  file_path <- paste0(input_dir, prop, "_200_hmm.txt")
  data <- read.table(file_path, header = FALSE, fill = TRUE)
  
  average_data <- colMeans(data)
  
  positions <- seq(-200, 200, length.out = ncol(data))
  
  average_df <- data.frame(Position = positions, Mean_Value = average_data)
  
  p <- ggplot(average_df, aes(x = Position, y = Mean_Value)) +
    geom_line(color = "blue") +
    theme_minimal() +
    labs(x = "Position", y = "Average HMM Feature Value", 
         title = paste("Average DNA Shape Prediction for controlled", prop, " in H.sapiens")) +
    scale_x_continuous(breaks = seq(-200, 200, 50))  
  plot_list_hsapiens_hmm[[prop]] <- p
}
grid.arrange(grobs = plot_list_hsapiens_hmm[1:5], ncol = 1)


human_opening_hmm <- read.table("~/Desktop/Promoters/dna_shape_control/hsapiens_Opening_200_hmm.txt", header = FALSE, fill = TRUE)
average_data <- colMeans(human_opening_hmm, na.rm = TRUE)
positions <- seq(-200, 200, length.out = 401)
average_df <- data.frame(Position = positions, Mean_Value = average_data)
ggplot(average_df, aes(x = Position, y = Mean_Value)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(x = "Position", y = "Average HMM Feature Value", 
       title = paste("Average DNA Shape Prediction for controlled", prop, " in H.sapiens")) +
  scale_x_continuous(breaks = seq(-200, 200, 50))

### Coefficient of variation in human data
human_buckle <- read.table("~/Desktop/Promoters/dna_shape/hsapiens_Opening_200.txt", header = FALSE)
human_buckle_hmm <- read.table("~/Desktop/Promoters/dna_shape_control/hsapiens_Buckle_200_hmm.txt", header = FALSE, fill = TRUE)

human_buckle <- na.omit(human_buckle)
human_buckle_hmm <- na.omit(human_buckle_hmm)

variation_coefficients <- apply(human_buckle, 2, function(x) sd(x) / mean(x))
filtered <- variation_coefficients[variation_coefficients < 200 & variation_coefficients > -200]

cv_df <- data.frame(
  position = 1:length(filtered), 
  variation_coefficient = filtered
)

ggplot(cv_df, aes(x = position, y = variation_coefficient)) +
  geom_point(color = "red", size = 1) +  
  labs(
    title = "Coefficient of Variation for Opening Parameter (H. sapiens)",
    x = "Position (Base Pairs)",
    y = "Coefficient of Variation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# Calculate control
variation_coefficients_control <- apply(human_buckle_hmm, 2, function(x) sd(x) / mean(x))
filtered_control <- variation_coefficients_control[variation_coefficients_control < 200 & variation_coefficients_control > -200]

cv_df_control <- data.frame(
  position = 1:length(filtered_control), 
  variation_coefficient = filtered_control
)

ggplot(cv_df_control, aes(x = position, y = variation_coefficient)) +
  geom_point(color = "red", size = 1) +  
  labs(
    title = "Coefficient of Variation for controlled Buckle Parameter (H. sapiens)",
    x = "Position (Base Pairs)",
    y = "Coefficient of Variation"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )


###one-way ANOVA 

human_buckle$promoter_id <- 1:nrow(human_buckle)
long_df <- melt(human_buckle, id.vars = "promoter_id", variable.name = "position", value.name = "value")

long_df$position <- as.numeric(gsub("V", "", long_df$position))

anova_results <- list()

for (i in 1:401) {
  pos_data <- subset(long_df, position == i)
  anova_test <- aov(value ~ factor(promoter_id), data = pos_data)
  anova_results[[i]] <- summary(anova_test)
  print(paste("Position", i, ": p-value =", summary(anova_test)[[1]][["Pr(>F)"]][1]))
}

p_values <- sapply(anova_results, function(x) x[[1]][["Pr(>F)"]][1])
pval_df <- data.frame(
  position = 1:401,
  p_value = p_values
)

ggplot(pval_df, aes(x = position, y = p_value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +  
  labs(title = "ANOVA p-values across positions", x = "Position", y = "p-value") +
  theme_minimal()
