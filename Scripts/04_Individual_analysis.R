##### Individual comparison #####

##### Paths #####
path_base <- normalizePath(getwd())
input_path <- normalizePath(file.path(path_base, "Input_Betas/"))
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

##### Load R packages #####
library(ggplot2)
library(tidyverse)
library(dplyr)
library(rstatix)
library(ggpubr)
library(patchwork)
library(ggbeeswarm)

##### Load results #####
load(file = file.path(plot_path, "data_clocks.RData"))
##### Function to get stats #####
get_stats <- function(df, variable_res) {
  formula_str <- paste(variable_res, "~ age_group")
  formula <- as.formula(formula_str)
  
  # Shapiro-Wilk test for normality by group
  p_shapiro <- df %>%
    group_by(age_group) %>%
    shapiro_test(!!sym(variable_res)) %>%
    pull(p) %>% min(na.rm = TRUE)
  
  if (p_shapiro < 0.05) {
    # Non-parametric test
    res_stats <- df %>% 
      wilcox_test(formula) %>%
      add_significance()
    
    message(paste(variable_res, ": Usando Wilcoxon Test"))
    
  } else {
    # Parametric test
    # Check Variances (Levene)
    p_levene <- df %>% levene_test(formula) %>% pull(p)
    
    if (p_levene > 0.05) {
      # Equal variances: Standard T-test
      res_stats <- df %>% 
        t_test(formula, var.equal = TRUE) %>%
        add_significance()
      message(paste(variable_res, ": Usando T-test (varianzas iguales)"))
    } else {
      # Unequal variances: Welch's T-test
      res_stats <- df %>% 
        t_test(formula, var.equal = FALSE) %>%
        add_significance()
      message(paste(variable_res, ": Usando Welch T-test (varianzas desiguales)"))
    }
  }
  
  # # Adjust position for plot
  res_stats <- res_stats %>% add_xy_position(x = "age_group", step.increase = 0.12)
  return(res_stats)
}


##### Style #####

#Boxplot
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
    legend.title = element_text(size = 18, face = "bold"), 
    legend.text = element_text(size = 16),                
    legend.key.size = unit(1.5, "cm"),                    
    )

# Lollyplot
estilo_publicacion2 <- estilo_publicacion + 
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    axis.minor.ticks.length.x = rel(0.4), 
    axis.ticks.length.x = unit(0.5, "cm"))
   
#Color palette
age_palette <- c(
  "≥ 12 years" = "#1E8E99",
  "< 12 years"  = "#FF8E32")

individual_palette <- c(
  "FA bone marrow" = "#3FA0FF", 
  "FA blood"       = "#D82632", 
  "FA age"  = "white"  
)

# Filter table for FA only
FA_results <- filter(df_resultados, Condition == "FA")

# Dataframe for segments 
df_segment <- FA_results %>%
  mutate(valor = age + EAAres_Horvath2) %>%
  group_by(Sample_Name) %>%
  summarise(
    age = dplyr::first(age),
    xend = max(c(age, valor), na.rm = TRUE),
    .groups = "drop"
  )
# Unique age per individual 
df_age_unique <- FA_results %>%
  group_by(Sample_Name) %>%
  summarise(age = dplyr::first(age))
# Plot
p_individual <- ggplot() +
# age group rectangles  
  annotate("rect", xmin = -Inf, xmax = 11.99,   
           ymin = -Inf, ymax = Inf,   
           fill = "#FF8E32", 
           alpha = 0.1) +              
  annotate("rect", xmin = 12, xmax = 28.5,   
           ymin = -Inf, ymax = Inf,   
           fill = "#1E8E99", 
           alpha = 0.1)+ 
  annotate("text", x = 5, y = 38, label = "Chronological age < 12 years", size = 5, fontface = "bold", angle = 0, color = "#FF8E32") +
  annotate("text", x = 20, y = 38, label = paste0("bold('Chronological age' >= 12 ~ 'years')"), size = 5, angle = 0, parse = TRUE, color = "#1E8E99")+
  
  # Segments up to maximun value
  geom_segment(
    data = df_segment,
    aes(
      x = 0,
      xend = xend,
      y = reorder(Sample_Name, -age),
      yend = reorder(Sample_Name, -age)
    ),
    color = "black",
    alpha = 0.6
  ) +
  # Chronological age
  geom_point(
    data = df_age_unique,
    aes(
      x = age,
      y = reorder(Sample_Name, -age),
      fill = "FA age"
    ),
    color = "black",
    stroke = 0.1,
    shape = 21,
    size = 4
  ) +
  # Chronological age + EAA in blood
  geom_point(
    data = FA_results %>% filter(Tissue == "blood"),
    aes(
      x = age + EAAres_Horvath2,
      y = reorder(Sample_Name, -age),
      fill = "FA blood"
    ),
    color = "black",
    stroke = 0.1,
    shape = 21,
    size = 4
  ) +
  
  # Chronological age + EAA in bone marrow
  geom_point(
    data = FA_results %>% filter(Tissue == "bone marrow"),
    aes(
      x = age + EAAres_Horvath2,
      y = reorder(Sample_Name, -age),
      fill = "FA bone marrow"
    ),
    color = "black",
    stroke = 0.1,
    shape = 21,
    size = 4
  ) +
  # Color setting
  scale_fill_manual(
    values = individual_palette,
    limits = c("FA age", "FA blood", "FA bone marrow"),
    labels = c("Chronological age", "CA + EAA in blood", "CA + EAA in bone marrow")
  ) +
  # Scales
  scale_x_continuous(breaks = seq(0, 50, by = 5), guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_discrete(expand = expansion(mult = c(0.025, 0.05)))+
  # Style
  estilo_publicacion2 +
  labs(
    title = "Horvath Skin & blood", 
    x = "Age in years", 
    y = "Patient ID",
    fill = "Age" 
  ) 

###### Age groups comparison (blood) #####

FA_blood <- filter(FA_results, Tissue == "blood")

FA_blood$age_group <- ifelse(FA_blood$age < 12, "< 12 years", "≥ 12 years")
# Define logical order
niveles_ordenados <- c("< 12 years", "≥ 12 years")
# Convert column to factor with specific order
FA_blood$age_group <- factor(FA_blood$age_group, 
                             levels = niveles_ordenados)

res_stats_blood <- get_stats(FA_blood, "EAAres_Horvath2")

p_blood <- ggplot(FA_blood, aes(x = age_group, y = EAAres_Horvath2, fill = age_group, show.legend = FALSE)) +
  scale_fill_manual(values = age_palette,
                    labels = c(expression("< 12 years"), 
                               expression(phantom() >= 12 ~ "years"))) +
  geom_boxplot(outlier.shape = NA, lwd = 0.5,median.linewidth= 0.5, width = 0.3, alpha = 0.2,  show.legend = FALSE) + 
  geom_quasirandom(aes(fill = age_group), size = 2,
                   shape= 21,  color = "black",stroke = 0.2,
                   width = 0.1, show.legend = FALSE) +
  scale_y_continuous(
    limits = c(NA, 30),                 
    breaks = seq(-5, 30, by = 5)      
  ) +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(
      x = stage(age_group, after_stat = as.numeric(x) + 0.3), 
      label = after_stat(round(y, 1))
    ), 
    vjust = 0.5,    
    hjust = 0,        
    color = "black", 
    fontface = "bold", 
    size = 4,
    inherit.aes = TRUE 
  ) +
  stat_pvalue_manual(res_stats_blood, label = "p", size = 5, tip.length = 0.01, fontface = "bold") +
  labs(title = "PB", x = NULL, y = "DNAm Age acceleration", fill = "") +
  estilo_publicacion +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.y = element_text(size = 16, face = "bold")) 

##### Age groups comparison (bone marrow) #####

FA_bm <- filter(FA_results, Tissue == "bone marrow")
FA_bm$age_group <- ifelse(FA_bm$age < 12, "< 12 years", "≥ 12 years")
# # Define logical order
niveles_ordenados <- c("< 12 years", "≥ 12 years")
# Transform column as factor
FA_bm$age_group <- factor(FA_bm$age_group, 
                          levels = niveles_ordenados)

res_stats_bm <- get_stats(FA_bm, "EAAres_Horvath2")

p_bm <- ggplot(FA_bm, aes(x = age_group, y = EAAres_Horvath2, fill = age_group)) +
  scale_fill_manual(
    values = age_palette,
    labels = c(
      expression(bold("< 12 years")), 
      expression(bold("" >= 12 ~ "years")) 
    )
  ) +
  geom_boxplot(outlier.shape = NA, lwd = 0.5,median.linewidth= 0.5, width = 0.3, alpha = 0.2) + 
  geom_quasirandom(aes(fill = age_group),
                   shape= 21,  color = "black",stroke = 0.2,
                   size = 2, width = 0.1, show.legend = FALSE) +
  scale_y_continuous(
    limits = c(NA, 30),                
    breaks = seq(0, 30, by = 5)      
  ) +
  stat_summary(
    fun = median, 
    geom = "text", 
    aes(
      x = stage(age_group, after_stat = as.numeric(x) + 0.3), 
      label = after_stat(round(y, 1))
    ), 
    vjust = 0.5,      
    hjust = 0,        
    color = "black", 
    fontface = "bold", 
    size = 4,
    inherit.aes = TRUE 
  ) +
  stat_pvalue_manual(res_stats_bm, label = "p", size = 5, tip.length = 0.01, fontface = "bold") +
  labs(title = "BM", x = NULL, y = "DNAm Age acceleration", fill = "") +
  estilo_publicacion +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 12, face = "bold"),  
    legend.title = element_text(size = 12, face = "bold"), 
    legend.key.size = unit(1, "cm"),
    axis.title.y = element_text(size = 16, face = "bold") 
    
  )
tissue_differences <- p_blood/p_bm+
  plot_layout(
    heights = c(4,4),
    guides = "keep",
    axis_titles = "collect")

Horvath2_analysis <- p_individual + inset_element(tissue_differences, left = 0.6, bottom = 0.2, right = 0.95, top = 0.95)  +
  plot_layout(
    guides = "keep")


# Save plot
ggsave(file.path(plot_path, "Horvath2_by_tissue_by_age.png"), 
       plot = Horvath2_analysis, 
       width = 15, 
       height = 9, 
       dpi = 300)