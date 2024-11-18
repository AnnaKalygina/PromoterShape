library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(viridis)
library(grid)
library(dplyr)
library(tidyr)


data <- read.csv("~/Desktop/Promoters/PromoterShape/combined_dna_shape_table.csv")
bootstrap_test <- read.csv("~/Desktop/Promoters/PromoterShape/bootstrap_experiment/raw_vs_hmm_bootstrap_test_cl_0.90_100_percent.csv")
data <- data[data$source %in% c("raw", "hmm", "probability_based"), ]

data[data$source == "probability_based", ]
bootstrap_test[bootstrap_test$significant == "True" & !is.na(bootstrap_test$upper_bound) & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185), ]

species <- c("athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae")
properties <-c( "MGW", "Buckle", "HelT", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt")
source <- c("raw","hmm", "probability_based")

pdf("~/Desktop/Promoters/PromoterShape/bootstrap_experiment/hmm_vs_raw_species_promoter_cl_0.90_100_percent.pdf", width = 8, height = 14)

for (spec in species) {
  for (prop in properties){
    df <- data[data$species == spec & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
    significance <- bootstrap_test[bootstrap_test$species == spec & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
    
    colors <- viridis(3)
    names(colors) <- c("raw", "hmm", "probability_based")
    
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
    
    significant_positions <- significance[significance$significant == "True", ]
    
    p3 <- ggplot(significance, aes(x = position, y = observed_difference)) +
      geom_line(size = 0.7) +  
      labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
      ylim(min(significance$raw_values), max(significance$raw_values)) +
      theme_minimal() +
      ggtitle("Difference between raw and hmm shapes") + 
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 14), 
        legend.position = "top" 
      ) +
      geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
                label = "*", color = "red", size = 5)
    
    
    combined_plot <- grid.arrange(p1, p2, p3, ncol = 1)
    print(combined_plot)
    #ggsave(filename = paste0("~/Desktop/Promoters/PromoterShape/raw_vs_control_graphs/", spec, "_", prop, ".pdf"), plot = combined_plot, width = 10, height = 8)
    
  }
}

dev.off()


### Build a bootstrap table
files <- list(
  "1" = "~/Desktop/Promoters/PromoterShape/bootstrap_experiment/raw_vs_hmm_bootstrap_test_cl_0.90_1_percent.csv",
  "10" = "~/Desktop/Promoters/PromoterShape/bootstrap_experiment/raw_vs_hmm_bootstrap_test_cl_0.90_10_percent.csv",
  "50" = "~/Desktop/Promoters/PromoterShape/bootstrap_experiment/raw_vs_hmm_bootstrap_test_cl_0.90_50_percent.csv",
  "100" = "~/Desktop/Promoters/PromoterShape/bootstrap_experiment/raw_vs_hmm_bootstrap_test_cl_0.90_100_percent.csv"
)


bootstrap_table <- data.frame()

for (sample_size in names(files)) {
  data <- read.csv(files[[sample_size]])
  print(sample_size)
  
  counts <- data %>%
    filter(significant == "True", position >= -185, position <= 185) %>%
    group_by(species, property) %>%
    summarise(significant_count = n(), .groups = "drop")
  
  counts <- counts %>%
    mutate(sample_size = as.numeric(sample_size))
  
  bootstrap_table <- rbind(bootstrap_table, counts)
}

complete_table <- bootstrap_table %>%
  complete(species, property, sample_size = c(1, 10, 50, 100), fill = list(significant_count = 0))

sorted_table <- complete_table %>%
  arrange(species, sample_size)

print(sorted_table)

write.csv(sorted_table, "~/Desktop/Promoters/PromoterShape/bootstrap_experiment/bootstrap_significance_summary.csv", row.names = FALSE)



cumulative_data <- bootstrap_table %>%
  group_by(species, sample_size) %>%
  summarise(cumulative_significant = sum(significant_count), .groups = "drop")

cumulative_significance <- ggplot(cumulative_data, aes(x = factor(sample_size), y = cumulative_significant, color = species)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(aes(group = species), linetype = "dotted") +
  scale_x_discrete(name = "Sample size, % from the total population", labels = c("1", "10", "50", "100")) +
  scale_y_continuous(name = "Significant areas count") +
  labs(title = "Cumulative significance count",
       color = "Species") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
png("~/Desktop/Promoters/PromoterShape/bootstrap_experiment/bootstrap_experiment_cumulative_significance.png", 
    width = 2000, height = 1200, res = 300)
print(cumulative_significance)
dev.off()

ten_percent_table <- complete_table[complete_table$sample_size == 10,]
fifty_percent_table <- complete_table[complete_table$sample_size == 50,]

ten <- ggplot(ten_percent_table, aes(x = property, y = significant_count, fill = property)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~species, scales = "free_y", ncol = 3) + 
  labs(
    title = "Significance count per property, 10% bootstrapping",
    x = "Property",
    y = "Significance count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_viridis_d()

fifty <- ggplot(fifty_percent_table, aes(x = property, y = significant_count, fill = property)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~species, scales = "free_y", ncol = 3) + 
  labs(
    title = "Significance count per property, 50% bootstrapping",
    x = "Property",
    y = "Significance count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_viridis_d()

combined_plot <- grid.arrange(ten, fifty, ncol = 2)
png("~/Desktop/Promoters/PromoterShape/bootstrap_experiment/bootstrap_property_wise_comparison.png", 
    width = 2000, height = 1000, res = 300)
print(combined_plot)
dev.off()

### Plotting shape exaggeration among species.
species <- c("athaliana", "celegans", "dmelanogaster", "hsapiens", "pfalciparum", "scerevisiae")
properties <-c( "MGW", "Buckle", "HelT", "Opening", "ProT", "Rise", "Roll", "Shear", "Shift", "Slide", "Stagger", "Stretch", "Tilt")
source <- c("raw","hmm")
pdf("shape_exaggeration_ss_10.pdf", width = 6, height = 18)

for (prop in properties){
  
  ### Athaliana
  df <- data[data$species == 'athaliana' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'athaliana' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  athaliana <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("A.thaliana shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  
  
  ### Dmelanogaster
  df <- data[data$species == 'dmelanogaster' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'dmelanogaster' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  dmelanogaster <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("D.melanogaster shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  
  ### Hsapiens
  df <- data[data$species == 'hsapiens' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'hsapiens' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  hsapiens <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("H.sapiens shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  ### Celegans
  df <- data[data$species == 'celegans' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'celegans' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  celegans <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("C.elegans shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  ### Scerevisiae
  df <- data[data$species == 'scerevisiae' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'scerevisiae' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  scerevisiae <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("S.cerevisiae shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  
  ### Pfalciparum
  df <- data[data$species == 'pfalciparum' & data$property == prop & (data$position >= -185 & data$position <= 185) , ]
  df_hmm <- df[df$source == "hmm",]
  df_raw <- df[df$source == "raw",]
  df_raw$difference <-df_raw$z_score - df_hmm$z_score
  
  significance <- bootstrap_test[bootstrap_test$species == 'pfalciparum' & bootstrap_test$property == prop & (bootstrap_test$position >= -185 & bootstrap_test$position <= 185) , ]
  significant_positions <- significance[significance$significant == "True", ]
  
  pfalciparum <- ggplot(significance, aes(x = position, y = observed_difference)) +
    geom_line(size = 0.7) +  
    labs(x = "Position around TSS (bp)", y = "Absolute difference value") +
    ylim(min(significance$raw_values), max(significance$raw_values)) +
    geom_vline(xintercept = significant_positions$position, size = 1, alpha = 0.1, color = 'red') +
    theme_minimal() +
    ggtitle(paste("P.falciparum shape exaggeration (raw - hmm) for", prop)) + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13), 
      axis.title = element_text(size = 13), 
      axis.text = element_text(size = 13), 
      legend.position = "top" 
    ) +
    geom_text(data = significant_positions, aes(x = position, y = observed_difference), 
              label = "*", color = "red", size = 5)
  
  combined_plot <- grid.arrange(athaliana, 
                                pfalciparum, 
                                scerevisiae, 
                                celegans, 
                                dmelanogaster, 
                                hsapiens, ncol = 1)
  print(combined_plot)
  
  
  
}

dev.off()

### Plotting all species table

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





