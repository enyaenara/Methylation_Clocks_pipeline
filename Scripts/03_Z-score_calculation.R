##### Script z-scores: Z-Score Calculation and Visualization #####

# This script calculates acceleration, standardizes it to z-score and generates a comparative 
# plot between patients (Fanconi) and controls to view the trend across clocks.

##### Load R packages #####
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
library(readxl)
library(patchwork)

##### Paths #####
path_base <- normalizePath(getwd())
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

##### Style and colors #####

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
    legend.key.size = unit(0.5, "cm"))                    

# Color palette
my_palette <- c(
  "FA bone marrow" = "#3FA0FF", 
  "FA blood"       = "#D82632",
  "control blood"  = "#40424f")

my_labels <- c(
  "control blood"  = "control PB",
  "FA blood"       = "FA PB",
  "FA bone marrow" = "FA BM"
)
##### Load necessary data #####
load(file = file.path(plot_path, "data_clocks.RData"))
RAFMex <- read_excel("Project_Data/Fanconi_anemia/RAFMex_Analysis.xlsx")

##### Standardization to Z-Score for EAA #####

# Create original vector with clock names
nombres_base <- c(
  "DNAmTL","DamAge",  
   "GrimAgeV2", "GrimAgeV1", 
  "PhenoAge","Hannum",
  "Horvath2", "Horvath1")

relojes_res <- paste0("EAAres_", nombres_base)

# Filter and calculate mean and SD of each clock 
# in control individuals
stats_controles <- df_resultados %>%
  filter(Group == "control blood") %>%
  summarise(across(all_of(relojes_res), 
                   list(mean = ~ mean(., na.rm = TRUE), 
                        sd = ~ sd(., na.rm = TRUE))))

# Save results in a new table
# where Z-scores will be stored

df_z_scores_res <- df_resultados

# Z-score calculation
for(reloj in relojes_res) {
  # Temporarily store mean and standard deviation of corresponding clock
  mean_ctrl <- stats_controles[[paste0(reloj, "_mean")]]
  sd_ctrl   <- stats_controles[[paste0(reloj, "_sd")]]
  
  # For each clock 
  # we apply formula Z-score = (x - control mean) / control sd
  df_z_scores_res[[paste0("z_", reloj)]] <- (df_z_scores_res[[reloj]] - mean_ctrl) / sd_ctrl
}

# Prepare data for plotting (Long Format)
# collapse Z-scores into three columns: Group, clock, and z-score
df_long_z_res <- df_z_scores_res %>%
  select(Group, starts_with("z_")) %>%
  pivot_longer(cols = starts_with("z_"), names_to = "Reloj", values_to = "Z_Score") %>%
  mutate(Reloj = gsub("z_EAAres_", "", Reloj))

# Calculate statistical summary by group and clock
df_summary_res <- df_long_z_res %>%
  group_by(Reloj, Group) %>%
  summarise(
    mean_z = mean(Z_Score, na.rm = TRUE),
    sd_z = sd(Z_Score, na.rm = TRUE),
    se_z = sd(Z_Score, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

# Change names for plotting 
df_summary_res <- df_summary_res %>%
  mutate(Reloj = recode_values(
    Reloj,
    "Horvath2" ~ "Skin & Blood Horvath",
    "Horvath1" ~ "Pan-Tissue Horvath",   
     default = Reloj                     
  ))

# To maintain correct order in plot
nombres_grafico <- c(
  "DNAmTL","DamAge",  
  "GrimAgeV2", "GrimAgeV1", 
  "PhenoAge","Hannum",
  "Skin & Blood Horvath", "Pan-Tissue Horvath")

# Force order with newly assigned names
df_summary_res$Reloj <- factor(df_summary_res$Reloj, levels = nombres_grafico)


# Plot
EAAres <- ggplot(df_summary_res, aes(x = Reloj, y = mean_z, color = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), 
                width = 0.2, position = position_dodge(0.5)) +
  geom_point(aes(fill = Group),position = position_dodge(0.5), size = 2, shape = 21, stroke = 0.1, color = "black") +
  coord_flip() + 
  estilo_publicacion +
  scale_color_manual(values = my_palette) +
  scale_fill_manual(values = my_palette) +
  labs(
    y = "Age acceleration (Z-Score units)",
    x = "Methylation Clock"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(plot_path,"EAAres_z-score_control.png"), plot = EAAres, width = 7, height = 6, dpi = 300)


###### Covariates Analysis #####

# RAFMex table preparation

RAFMex$Sample_Name <- gsub("-", "", RAFMex$Sample_ID) 
RAFMex$tissue <- sub("whole blood", "blood", RAFMex$tissue)

# Since both tables have the same Sample_Name, we join them
resultados_covariables <- left_join(df_resultados, RAFMex, by = "Sample_Name")
# Make it correctly detect NAs
resultados_covariables$Genotype <- ifelse(resultados_covariables$Genotype == "NA",NA, resultados_covariables$Genotype)
# For genotype, set non-FANCA to "Other" category
resultados_covariables$Genotype <- case_when(
  is.na(resultados_covariables$Genotype) ~ NA_character_,
  resultados_covariables$Genotype == "FANCA" ~ "FANCA",
  TRUE ~ "Other" 
)
# Assign folder where results will be saved
path_covariables <- file.path(plot_path, "Analisis_por_Covariable")

# Define column names to covariates 
covariables <- c("sex.x", "Transfusion_dependency", "Genotype", "Androgen_treatment")

# Adjust names for plotting
nombres_limpios <- c(
  "sex.x"                  = "Sex",
  "Genotype"               = "Genotype",
  "Transfusion_dependency" = "Transfusion Dependency",
  "Androgen_treatment"     = "Androgen treatment"
)

# Loop to generate plots for each covariate
for (covar in covariables) {
  
  message(paste("Ajustando por covariable:", covar))
  # Keep all controls
  # and patients having a valid value in the covariate
  df_temp <- resultados_covariables %>%
    filter(Group == "control blood" | !is.na(.data[[covar]]))
  
  # Sub-loop for Clock adjustment
  for(reloj in relojes_res) {
    
    df_controles <- df_temp %>% filter(Group == "control blood")
    rownames_controles <- rownames(df_controles)
    # Create vector with original EAAres_ values
    EAA_final <- df_temp[[reloj]]
    names(EAA_final) <- rownames(df_temp)
    
    if (covar == "sex.x") {
      # Sex: Varies in controls and patients. 
      # Controls are taken as reference
      modelo_ajuste <- lm(as.formula(paste(reloj, "~ sex.x")), data = df_controles)
      # Adjustment (Observed - Predicted)
      valor_esperado <- predict(modelo_ajuste, newdata = df_temp)
      EAA_final <- df_temp[[reloj]] - valor_esperado
      
    } else {
      # Clinical variables: Vary only in patients.
      df_pacientes <- df_temp %>% filter(Group != "control blood")
      
      if (n_distinct(df_pacientes[[covar]], na.rm = TRUE) < 2) {
        stop(paste("La covariable", covar, "no tiene suficientes niveles. Revisar datos."))
      }
      
      # Train model to see covariate effect in patients
      modelo_ajuste <- lm(as.formula(paste(reloj, "~", covar)), data = df_pacientes)
      
      # Pure clinical effect is total prediction minus basal intercept
      efecto_clinico <- predict(modelo_ajuste) - coef(modelo_ajuste)[1]
      
      # Subtract excess years caused by covariate from patients.
      # Controls remain intact (clinical effect is 0 as they don't have the variable)
     EAA_final[rownames(df_pacientes)] <- df_pacientes[[reloj]] - efecto_clinico
    }
    
    ##### Z-Score calculation with controls as reference #####
    # Take control reference values
    valores_referencia_control <- EAA_final[df_temp$Group == "control blood"]
    
    # Define control mean and sd
    mean_control <- mean(valores_referencia_control, na.rm = TRUE)
    sd_control   <- sd(valores_referencia_control, na.rm = TRUE)
    
    # Final Z-score
     z_EAA_adj <- (EAA_final - mean_control) / sd_control
    
    df_temp[[paste0("z_adj_", reloj)]] <- z_EAA_adj
  } # End of sub-loop
  
  ##### Prepare data for plotting #####
  # Convert to long format
  df_long_temp <- df_temp %>%
    select(Group, starts_with("z_adj_")) %>%
    pivot_longer(cols = starts_with("z_adj_"), names_to = "Reloj", values_to = "Z_Adj") %>%
    mutate(Reloj = gsub("z_adj_EAAres_", "", Reloj))
  
  # Get mean and error by group
  df_summary_temp <- df_long_temp %>%
    group_by(Reloj, Group) %>%
    summarise(
      mean_z = mean(Z_Adj, na.rm = TRUE),
      se_z = sd(Z_Adj, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Adjust names for plotting
  df_summary_temp <- df_summary_temp %>%
    mutate(Reloj = case_match(
      Reloj,
      "Horvath2" ~ "Skin & Blood Horvath",
      "Horvath1" ~ "Pan-Tissue Horvath",   
      .default = Reloj                     
    ))
  
  # To maintain correct order in plot
  nombres_grafico <- c(
    "DNAmTL","DamAge",  
    "GrimAgeV2", "GrimAgeV1", 
    "PhenoAge","Hannum",
    "Skin & Blood Horvath", "Pan-Tissue Horvath")
  
  # Force order with newly assigned names
  df_summary_temp$Reloj <- factor(df_summary_temp$Reloj, levels = nombres_grafico)
  
  ### Adjust annotations in plots ###
  
  # Define rightmost point (max Z + its error bar)
  limite_derecho <- max(df_summary_temp$mean_z + df_summary_temp$se_z, na.rm = TRUE)
  
  if (covar == "sex.x") {
    posicion_anotacion_y <- limite_derecho + 0.5
    posicion_texto_y     <- posicion_anotacion_y + (limite_derecho * 0.05)}
  
  else {posicion_anotacion_y <- limite_derecho + 0.5
  posicion_texto_y     <- posicion_anotacion_y + 0.2}
  
  ##### Plot by covariate #####
  # dotted line represents controls
 p <- ggplot(df_summary_temp%>% filter(Group != "control blood"), aes(x = Reloj, y = mean_z, color = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_errorbar(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), 
                  width = 0.2, position = position_dodge(0.5)) +
    geom_point(aes(fill = Group), position = position_dodge(0.5), 
               size = 2, shape = 21, stroke = 0.1, color = "black") +
    coord_flip() +
    estilo_publicacion +
    scale_color_manual(values = my_palette, labels = my_labels) +
    scale_fill_manual(values = my_palette , labels = my_labels) +
    labs(
      title = paste("Covariate: ", nombres_limpios[covar]),
      y = "Age acceleration (z-score units)",
      x = "Methylation Clock"
    ) +
    theme(legend.position = "right") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    # Annotations for first generation clocks
    annotate("segment", x = 5.7, xend = 8.5, y = posicion_anotacion_y, yend = posicion_anotacion_y, size = 1.2, color = "black") + 
    annotate("text", x = 7, y = posicion_texto_y, label = "1st Generation", size = 3, fontface = "bold", angle = -90) +
    # Annotations for second generation clocks
    annotate("segment", x = 0.7, xend = 5.3, y = posicion_anotacion_y, yend = posicion_anotacion_y, size = 1.2, color = "black") + 
    annotate("text", x = 3, y = posicion_texto_y, label = "2nd Generation", size = 3, fontface = "bold", angle = -90) 
  
  # Assign object with unique name to R environment
  nombre_objeto_r <- paste0("plot_", covar)
  assign(nombre_objeto_r, p)
}
