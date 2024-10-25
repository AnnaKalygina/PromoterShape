library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(viridis)


data <- read.csv("~/Downloads/combined_dna_shape_table.csv")

species <- c("athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae")
properties <- c("Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt")
sources <- c("raw", "shuffled", "hmm")


for (spec in species) {
  for (prop in properties){
    df <- data[data$species == spec & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
    
    colors <- viridis(3)
    names(colors) <- c("raw", "hmm", "shuffled")
    
    p1 <- ggplot(df, aes(x = position, y = value)) +
      geom_point(aes(color = source), size = 1, alpha = 0.7) +  
      scale_color_manual(values = colors) + 
      labs(x = "Position around TSS (bp)", y = "Average Z-score") +
      facet_wrap(~source, ncol = 3) + 
      ggtitle(paste(prop, "Z-score for", spec)) +  
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "none",
        panel.spacing = unit(22, "pt")
      )
    
    p2 <- ggplot(df, aes(x = position, y = value, color = source)) +
      geom_line()+  
      scale_color_manual(values = colors) + 
      labs(x = "Position around TSS (bp)", y = "Average Z-score") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "top"  # Show the legend in overlay plot
      )
    
    combined_plot <- grid.arrange(p1, p2, ncol = 1)
    ggsave(filename = paste0("~/Downloads/", spec, "_", prop, ".pdf"), plot = combined_plot, width = 10, height = 8)
    
  }
}