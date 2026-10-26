# Analysis of comparison with progeroids


##### Load R packages #####
library(methylCIPHER) #HigginsChenLab/methylCIPHER
library(limma)
library(ggplot2)
library(ggplotify)
library(ggpubr)
library(minfi)
library(readr)
library(tidyverse)
library(dplyr)
library(rstatix)
library(patchwork)
library(EpiDISH)
library(pheatmap)
library(cowplot)

##### Paths #####
path_base <- normalizePath(getwd())
input_path <- normalizePath(file.path(path_base, "Input_Betas"))
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

##### Style #####
estilo_publicacion <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    axis.ticks.length = unit(0.2, "cm"),                        
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10),                
    legend.key.size = unit(1, "cm"),                    
  )

# Color palette
progeroid_palette <- c(
  "control blood"  = "#40424f", 
  "FA blood"       = "#D82632", 
  "FA bone marrow" = "#3FA0FF", 
  "WS blood"       ="#E69F00",
  "HGPS blood"     ="#009E73",
  "NCPL blood"     ="#CC79A7"
)

# Plot labels
progeroid_labels <- c(
  "control blood"  = "control PB",
  "FA blood"       = "FA PB",
  "FA bone marrow" = "FA BM",
  "WS blood"       ="WS PB",
  "HGPS blood"     ="HGPS PB",
  "NCPL blood"     ="NCPL PB"
)

##### Dataframe preparation #####
# Load results obtained from processing
beta_harmonized <- read_csv(file.path(input_path,"beta_harmonized_completo.csv"))
Pheno <- read_csv(file.path(input_path,"all_targets_completo.csv"))

# Convert to dataframe and set first column as row names
beta_harmonized <- as.data.frame(beta_harmonized) %>%
  column_to_rownames(var = "...1")

#Adjust pheno
Pheno <- as.data.frame(Pheno) 
Pheno$Tissue <- sub("whole blood", "blood", Pheno$Tissue)
Pheno$Group <- paste(Pheno$Condition, Pheno$Tissue, sep = " ")

# Adjust basename to match names in df_final
Pheno$Basename <- basename(Pheno$Basename)

# Adjust name only for external controls
Pheno$Sample_ID <- ifelse(
  grepl("GSM", Pheno$Sample_Name) & 
    Pheno$Condition == "control" & 
    # Controls from these datasets (WS and HGPS) don't need this adjustment
    # since we saved complete targets
    !grepl("^GSM(269|55)", Pheno$Sample_Name), 
  paste(Pheno$Basename, "control"), 
  Pheno$Basename
)

#Make sure betas match Pheno
df_final <- select(beta_harmonized, Pheno$Sample_ID)

##### Epigenetic clocks analysis #####
# Transpose data frame to be compatible with methylCIPHER
df_final <- t(df_final)

# Convert variables to factors
Pheno$Group <- as.factor(Pheno$Group)
Pheno$Tissue    <- as.factor(Pheno$Tissue)


# Check missing sites
missing_CpG <- getClockProbes(df_final)
missing_CpG$Percent.Present <- sub("%","", missing_CpG$Percent.Present)
missing_CpG$Percent.Present <- as.numeric(missing_CpG$Percent.Present)

# Run specifically for GrimAge
Pheno$Female <- ifelse(Pheno$sex == "F", 1, 0)
Pheno$Age <- Pheno$age

# Define clocks to calculate
userClocks <- c("calcHorvath1", "calcHorvath2","calcHannum",
                "calcDNAmTL", "calcDamAge","calcPhenoAge"
)

# Calculate all clocks at once (except GrimAge)
df_resultados <- calcUserClocks(
  userClocks,
  df_final,
  Pheno,
  imputation = FALSE
)
# Calculate both versions of GrimAge
GrimAgeV1 <- calcGrimAgeV1(df_final, Pheno, ID = "Sample_ID")
GrimAgeV2 <- calcGrimAgeV2(df_final, Pheno, ID = "Sample_ID")

# Merge with former results
df_resultados$GrimAgeV1 <- GrimAgeV1$GrimAgeV1
df_resultados$GrimAgeV2 <- GrimAgeV2$GrimAgeV2

# Get clock names (without 'calc')
relojes <- c(gsub("calc", "", userClocks), "GrimAgeV1", "GrimAgeV2")


##### Analysis #####

# Loop to generate regression plots

# Order table
df_resultados_ordenado <- df_resultados %>%
  arrange(Group)

relojes <- c("Horvath1","Hannum","PhenoAge",
            "GrimAgeV1","GrimAgeV2","DamAge", "DNAmTL" 
)

lista_de_plots <- list()
lista_de_stats <- list()

for (reloj in relojes) {
  cat("Procesando:", reloj, "\n")
  # Adjust title name
  name <- case_when(
    reloj == "Horvath1" ~ "Horvath Pan-tissue",
    TRUE                       ~ reloj 
  )
  
  # Adjust Y-axis name 
  eje_y_label <- if (reloj == "DNAmTL") {
    "Predicted Telomere Lenght"
  } else {
    "DNAm Age" 
  }
 
  # Plot
  try({
    p <- ggplot(df_resultados_ordenado, aes(x = age, y = .data[[reloj]], color = Group)) +
      scale_color_manual(values = progeroid_palette, labels = \(x) progeroid_labels[x])+ 
      geom_point(aes(fill = Group),
                 color = "black",
                 stroke = 0.2,
                 shape = 21,
                 size = 3) +
      scale_fill_manual(values = progeroid_palette, labels = \(x) progeroid_labels[x])+
      geom_smooth(aes(group = Group), method = "lm", linewidth = 0.8, se = FALSE) +
      labs(title = name, x = "Chronological Age",
           y = eje_y_label,
           fill = "Group", 
      ) +
      estilo_publicacion +
      theme(
        axis.title.y = element_text(size=16),
        legend.position = "bottom")
    
    lista_de_plots[[reloj]] <- p
    
    ggsave(file.path(plot_path,"Analisis_vs_progeroides", paste0("Plot_", reloj, ".png")), plot = p,
           width = 7,
           height = 6,
           dpi = 300)
  }, silent = FALSE)
}

cat("Scatter plots completos. Revisa la carpeta 'Plot_clocks/Analisis_vs_progeroides'.\n")

p_leyenda <- wrap_plots(lista_de_plots, nrow = 4) + 
  plot_layout(guides = "collect") & theme(legend.position = "right")
leyenda_a <- cowplot::get_legend(p_leyenda)

lista_sin_leyenda <- lapply(lista_de_plots, function(p) p + theme(legend.position = "none"))

# Add legend
lista_sin_leyenda[[8]] <- wrap_elements(leyenda_a)
p_relojes <- wrap_plots(lista_sin_leyenda, nrow = 4) 
 
p_a <- wrap_elements(p_relojes) + labs(tag = "A") &
  theme(plot.tag = element_text(size = 18, face = "bold"))
  
##### Panel B: Heatmap Episignature #####

# Get from Script 08_FA_WS_HGPSvsControls
annotation_col <- data.frame(
  Group = factor(df_plot_firma$Group))

rownames(annotation_col) <- colnames(beta_firma_user)

ann_colors <- list(
  Group = progeroid_palette
)

hm_obj <- pheatmap(beta_firma_user,
                   clustering_method = "ward.D2",
                   annotation_col = annotation_col,    
                   annotation_colors = ann_colors,     
                   annotation_legend = TRUE,
                   show_rownames = FALSE,
                   show_colnames = FALSE,              
                   fontsize = 12,
                   color = colorRampPalette(c("blue", "white", "yellow"))(100),
                   silent = TRUE)

hm_signature <- as.ggplot(hm_obj)
  
hm_signature <- wrap_elements(hm_signature) + labs(tag = "B")  &
  theme(plot.tag = element_text(size = 18, face = "bold"))

p_b <- plot_spacer()/hm_signature/plot_spacer() + plot_layout(heights = c(4,10,4))

##### Multidimensional analysis #####
p_c <- ggplot(df_plot_firma, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(fill = Group), color = "black",
             stroke = 0.2,
             shape = 21, size = 3) +
  scale_fill_manual(values = progeroid_palette) +
  scale_y_continuous(expand = expansion(mult = 0.15)) + 
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  labs(x = "MDS 1", y = "MDS 2") +
  estilo_publicacion

p_c <- p_c + labs(tag = "C")

episignature <- p_b/p_c + plot_layout(nrow = 2)
progeroid_analysis <- wrap_plots(p_a,episignature, ncol=2, widths = c(9,15))


# Save
ggsave(file.path(plot_path, "Progeroid_analysis.pdf"), 
       plot = progeroid_analysis, 
       width = 24,   
       height = 18,  
       units = "cm",
       dpi = 300,
       scale= 2,
       useDingbats = FALSE)


cat("Guardar datos")
save(missing_CpG, df_final, Pheno, df_resultados, df_plot_firma, beta_firma_user,
        file = file.path(plot_path, "Analisis_vs_progeroides/data_clocks_progeria.RData"), compress = TRUE)

