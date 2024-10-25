library(ggplot2)
library(ggpubr)


input_dir <- "~/Downloads/"
species <- c("athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae")
buckle <- data.frame(matrix(ncol = length(species) + 1, nrow = 401))
colnames(buckle) <- append("position", species)
buckle$position <- seq(-200,200, length.out = 401)

for (spec in species) {
  file_path <- paste0(input_dir, spec, "_Buckle_200.txt")
  data <- read.table(file_path, header = FALSE)
  

  data_z_scored_rows <- t(apply(data, 1, function(x) (x - mean(x)) / sd(x)))
  averaged_data <- colMeans(data_z_scored_rows)
  buckle[[spec]] <- averaged_data
  

  plot_data <- data.frame(position = buckle$position[26:376], z_score = buckle[[spec]][26:376])
  pl <- ggplot(plot_data, aes(x = position, y = z_score)) +
    geom_line(color = "blue") +  
    labs(x = "Position around TSS (bp)", y = "Average Z-score") +
    ggtitle(paste("Normalised Buckle in", spec)) + 
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  assign(paste0(spec, "_plot"), pl)
}


ggarrange(get("athaliana_plot"), get("hsapiens_plot"), 
          get("dmelanogaster_plot"), get("celegans_plot"), 
          get("pfalciparum_plot"), get("scerevisiae_plot"),
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)


### Plotting pairwise comparison

par(mfrow=c(1, 1))
matplot(buckle$position[26:376], cbind(buckle$hsapiens[26:376], buckle$dmelanogaster[26:376]), 
        type = "l", lty = 1, lwd = 2,
        col = c("blue", "green"), xlab = "Position around TSS (bp)", 
        ylab = "Average z-score", main = "Comparison of buckle prediction",
        cex.lab = 1.2,  
        cex.axis = 1.1, 
        cex.main = 1.5)
legend("bottomright", legend = c("H. sapiens", "D. melanogaster"), 
       col = c("blue", "green"), 
       lty = 1,
       cex = 0.8)


par(mfrow = c(1, 1))

# Plot the buckle predictions with custom line thickness and colors
matplot(buckle$position[26:376], 
        cbind(buckle$hsapiens[26:376], buckle$dmelanogaster[26:376]), 
        type = "l", lty = 1, lwd = 2,  # Line thickness
        col = c("blue", "green"), 
        xlab = "Position around TSS (bp)", 
        ylab = "Average z-score", 
        main = "Comparison of buckle prediction",
        cex.lab = 1.2,  # Increase axis label size
        cex.axis = 1.1,  # Increase axis tick size
        cex.main = 1.5)  # Increase title size

# Add a smaller legend with custom size
legend("topright", legend = c("H. sapiens", "D. melanogaster"), 
       col = c("blue", "green"), 
       lty = 1, lwd = 2,  # Line thickness in legend
       cex = 0.8)  # Reduce the size of the legend
