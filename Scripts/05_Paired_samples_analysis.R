##### Comparison between FA with both tissues from the same sampling age #####


##### Load R packages #####
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)
library(ggplotify)
library(rstatix)
library(tidyr)
library(viridis)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(readxl)
library(pheatmap)
library(ggbeeswarm)

##### Paths #####
path_base <- normalizePath(getwd())
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

###### Load results #####
load(file = file.path(plot_path, "ALL_FA_vs_ALL/data_clocks.RData"))

##### Style #####
estilo_publicacion <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),                        
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 10, face = "bold"), 
    legend.text = element_text(size = 8),               
    legend.key.size = unit(0.5, "cm"),                    
    )

##### Color palette #####
tissue_palette <- c(
  "bone marrow" = "#3FA0FF",
  "blood"       = "#D82632")
##### Function to perform paired analysis #####
get_paired_stats <- function(df, variable_res) {
  formula_str <- paste(variable_res, "~ Tissue")
  formula <- as.formula(formula_str)
  # Normality is evaluated on the difference between pairs
  diferencias <- df %>%
    select(Sample_Name, Tissue, !!sym(variable_res)) %>%
    tidyr::pivot_wider(names_from = Tissue, values_from = !!sym(variable_res)) %>%
    mutate(diff = .[[2]] - .[[3]]) 
  
  p_shapiro <- shapiro_test(diferencias, diff) %>% pull(p)
  
  if (p_shapiro < 0.05) {
    # Non-parametric test: Wilcoxon Signed-Rank Test
    result <- df %>% 
      wilcox_test(formula, paired = TRUE) %>% 
      add_significance() 
    message(paste(variable_res, ": Usando Wilcoxon signed-rank test"))
    
  } else {
    # Parametric test: Paired t-test
    result <- df %>% 
      t_test(formula, paired = TRUE) %>% 
      add_significance()
    message(paste(variable_res, ": Usando Paired t-test"))
  }
  
  # Adjust position
  result <- result %>% add_xy_position(x = "Tissue", step.increase = 0.12)
  return(result)
}

##### Dataframe preparation #####

# Filter dataframe for FA patients 
FA_results <- filter(df_resultados, Condition == "FA")
# Filter those with both tissues and same samplig age
FA_filtered <- FA_results %>%
  group_by(Sample_Name) %>%
  filter(n() == 2) %>%
  group_by(Sample_Name) %>%
  filter(n_distinct(age) == 1) %>%
  ungroup()
# Vector with results to analyze
relojes <- c("EAAres_Horvath1", "EAAres_Horvath2", 
                "EAAres_Hannum",  "EAAres_DNAmTL", 
               "EAAres_DamAge","EAAres_PhenoAge", 
             "EAAres_GrimAgeV1", "EAAres_GrimAgeV2")

lista_de_plots <- list()

for (reloj in relojes) {
  cat("Calculando estadísticos para:", reloj, "\n")
  
  try({
    res_stats <- get_paired_stats(FA_filtered, reloj)
    
    # Adjust title names
    name <- case_when(
      reloj == "EAAres_Horvath1" ~ "Horvath Pan-tissue",
      reloj == "EAAres_Horvath2" ~ "Horvath Skin & blood",
      TRUE                       ~ sub("EAAres_", "", reloj)
    )
    
    # Adjust Y axis labels
    eje_y_label <- if (reloj == "EAAres_DNAmTL") {
      "DNAmTL residual"
    } else {
      "DNAm age acceleration" 
    }
    
    boxplot <- ggplot(FA_filtered, aes(x = Tissue, y = .data[[reloj]])) +
      scale_fill_manual(values = tissue_palette, labels = c("Periphereal blood", "Bone marrow")) +
      geom_boxplot(aes(fill = Tissue), outlier.shape = NA, lwd = 0.2, 
                   median.linewidth = 0.5, width = 0.4, alpha = 0.2, show.legend = TRUE) + 
      geom_line(aes(group = Sample_Name), color = "gray", linewidth = 0.3) +
      geom_point(aes(fill = Tissue), size = 1.2, shape = 21, color = "black", 
                 stroke = 0.2, show.legend = FALSE) +
      stat_summary(
        fun = median, 
        geom = "text", 
        aes(
          x = stage(Tissue, after_stat = as.numeric(x) + 0.3), 
          label = after_stat(round(y, 1))
        ), 
        vjust = 0.5,      
        hjust = 0,        
        color = "black", 
        fontface = "bold", 
        size = 2, 
        inherit.aes = TRUE 
      ) +
      stat_pvalue_manual(res_stats, label = "p", size = 2.5, tip.length = 0.05) +
      labs(title = name, y = eje_y_label, x = NULL, fill = "Tissue") + 
      estilo_publicacion +
      theme(legend.position = "bottom",
            legend.key.size = unit(0.3, "cm"),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    
    # Save to plot list 
    lista_de_plots[[reloj]] <- boxplot
    
  }, silent = FALSE)
}


# Wrap plots
plot_final <- wrap_plots(lista_de_plots, ncol = 4) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save plot
ggsave(
  file.path(plot_path, "Paired_tissue_analysis.png"), 
  plot = plot_final, 
  width = 6, 
  height = 4, 
  dpi = 300,
  scale = 2
)