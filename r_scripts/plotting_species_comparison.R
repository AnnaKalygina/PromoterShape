library(ggplot2)
library(dplyr)
library(gridExtra)


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

data <- read.csv("../tennisnyjmac/Desktop/Promoters/PromoterShape/combined_dna_shape_table.csv")
cosine_similarity_data <- read.csv("~/Desktop/Promoters/PromoterShape/cosine_shapes_similarity.csv")

filtered_data <- data %>% 
  filter(source == "raw", position >= -185, position <= 185)

species_list <- unique(filtered_data$species)
property_list <- unique(filtered_data$property)

pdf("Species_Shapes_Comparison_Plots.pdf", width = 14, height = 8)


for (prop in property_list) {
  prop_data <- filtered_data %>% filter(property == prop)
  plot_list <- list()

  species_pairs <- combn(species_list, 2, simplify = FALSE)
  
  for (pair in species_pairs) {
    species1 <- pair[1]
    species2 <- pair[2]
    
    species1_data <- prop_data %>% filter(species == species1)
    species2_data <- prop_data %>% filter(species == species2)
    
    comparison_data <- data.frame(
      position = species1_data$position,
      value1 = species1_data$value,
      value2 = species2_data$value
    )
    
    cosine_similarity_value <- cosine_similarity_data %>%
      filter(species1 == !!species1, species2 == !!species2, property == prop) %>%
      pull(cosine_similarity)
    
    plot <- ggplot(comparison_data, aes(x = position)) +
      geom_line(aes(y = value1, color = species1), size = 0.6) +
      geom_line(aes(y = value2, color = species2), size = 0.6) +
      labs(title = paste("Comparison of", prop, "between", species1, "and", species2),
           x = "Position around TSS (bp)",
           y = "Average z-score") +
      scale_color_manual(values = c("blue", "green"), labels = c(species1, species2)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 7),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 6)) +
      annotate("text", x = 180, y = 0.25,
               label = round(cosine_similarity_value, 2),
               color = "black", fontface = "bold", hjust = 1.6, size = 4)
    
    plot_list[[paste(species1, species2)]] <- plot
  }
  
  grid_plot <- marrangeGrob(plot_list, nrow = 2, ncol = 2, top = paste("Property:", prop))
  print(grid_plot)
}

dev.off()


