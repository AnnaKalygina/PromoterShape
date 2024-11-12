library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(viridis)
library(grid)


data <- read.csv("~/Desktop/Promoters/PromoterShape/combined_dna_shape_table.csv")

species <- c("athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae")
properties <-c( "MGW", "Buckle", "HelT", "MGW", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt")
sources <- c("raw", "shuffled", "hmm")

pdf("hmm_vs_raw_species_prmoter.pdf", width = 8, height = 14)

for (spec in species) {
  for (prop in properties){
    df <- data[data$species == spec & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
    
    colors <- viridis(3)
    names(colors) <- c("raw", "hmm", "shuffled")
    
    p1 <- ggplot(df, aes(x = position, y = z_score)) +
      geom_point(aes(color = source), size = 1, alpha = 0.7) +  
      scale_color_manual(values = colors) + 
      labs(x = "Position around TSS (bp)", y = "Average Z-score") +
      facet_wrap(~source, ncol = 3) + 
      ggtitle(paste(prop, "Z-score for", spec)) +  
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "none",
        panel.spacing = unit(22, "pt")
      )
    
    p2 <- ggplot(df, aes(x = position, y = z_score, color = source)) +
      geom_line()+  
      scale_color_manual(values = colors) + 
      labs(x = "Position around TSS (bp)", y = "Average Z-score") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "top" 
      )
    
    df_hmm <- df[df$source == "hmm",]
    df_raw <- df[df$source == "raw",]
    df_raw$difference <-df_raw$z_score - df_hmm$z_score
    
    p3 <- ggplot(df_raw, aes(x = position, y = difference)) +
      geom_line(size = 0.7)+  
      labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
      ylim(min(df_raw$z_score), max(df_raw$z_score)) +
      theme_minimal() +
      ggtitle("Difference between raw and hmm shapes") + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "top" 
      )
    
    combined_plot <- grid.arrange(p1, p2, p3, ncol = 1)
    print(combined_plot)
    #ggsave(filename = paste0("~/Desktop/Promoters/PromoterShape/raw_vs_control_graphs/", spec, "_", prop, ".pdf"), plot = combined_plot, width = 10, height = 8)
    
  }
}

dev.off()

### Plotting all properties in one file

species_list <- unique(data$species)
shape_features <- unique(data$property)

for (spec in species_list) {
  
  species_data <- data[data$species == spec & data$source == "raw" & (data$position >= -185 & data$position <= 185), ]
  
  plots <- list()
  
  for (feature in shape_features) {
    print(feature)
    plot <- ggplot(species_data[species_data$property == feature,],
                   aes(x = position, y = value)) +
      geom_line(color = "blue") + 
      labs(title = feature,
           x = "Position around TSS (bp)",
           y = "Mean Z-score value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.title = element_text(),
            axis.text = element_text(size = 7),
            panel.border = element_rect(color = "black", fill = NA, size = 0.7))
    
    plots[[feature]] <- plot
  }
  
  
  grid_plot <- arrangeGrob(grobs = plots, ncol = 4,
                           top = textGrob(paste0("Shape features graphs for ", spec),
                                          gp = gpar(fontsize = 20, fontface = "bold")))
  
  #ggsave(paste0(spec, "_shape_plots.pdf"), grid_plot, width = 16, height = 14)
}



