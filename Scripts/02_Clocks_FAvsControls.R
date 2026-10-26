##### Load R packages #####
library(methylCIPHER) #HigginsChenLab/methylCIPHER
library(limma)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(minfi)
library(readr)
library(tidyverse)
library(dplyr)
library(rstatix)
library(patchwork)

##### Paths #####
path_base <- normalizePath(getwd())
input_path <- normalizePath(file.path(path_base, "Input_Betas"))
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

##### Style configuration #####
# We defined this style for the paper figures, thus it is optional
estilo_publicacion <- theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    axis.line = element_line(color = "black", linewidth = 0.8), 
    axis.ticks = element_line(color = "black", linewidth = 0.8), 
    axis.ticks.length = unit(0.2, "cm"),                        
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10),                
    legend.key.size = unit(2, "cm"),                    
    )

# Color palette for group comparisons
my_palette <- c(
  "control blood"  = "#40424f", 
  "FA blood"       = "#D82632", 
  "FA bone marrow" = "#3FA0FF"  
)

# Labels appearing in the plots
my_labels <- c(
  "control blood"  = "control PB",
  "FA blood"       = "FA PB",
  "FA bone marrow" = "FA BM"
)

##### Input files preparation #####

# Load results obtained from processing
beta_harmonized <- read_csv(file.path(input_path,"beta_harmonized.csv"))
Pheno <- read_csv(file.path(input_path,"all_targets.csv"))

# Convert to data frame and set first column as row names
beta_harmonized <- as.data.frame(beta_harmonized) %>%
  column_to_rownames(var = "...1")

# Adjust Pheno so variables are easy to plot
Pheno <- as.data.frame(Pheno) 
Pheno$Group <- paste(Pheno$Condition, Pheno$Tissue, sep = " ")

# Modify basename to match names in df_final
Pheno$Basename <- basename(Pheno$Basename)
Pheno$Sample_ID <- ifelse(grepl("GSM", Pheno$Sample_Name), 
                         paste(Pheno$Basename, "control"), 
                         Pheno$Basename)
# Filter betas to ensure match with Pheno
df_final <- select(beta_harmonized, Pheno$Sample_ID)

# Transpose data frame to be compatible with methylCIPHER
df_final <- t(df_final)

# Convert variables to factors
Pheno$Group <- as.factor(Pheno$Group)
Pheno$Tissue    <- as.factor(Pheno$Tissue)

##### Calculation of epigenetic clock ages for both tissues #####

# Verify that at least 80% of CpG sites are present
missing_CpG <- getClockProbes(df_final)
missing_CpG$Percent.Present <- sub("%","", missing_CpG$Percent.Present)
missing_CpG$Percent.Present <- as.numeric(missing_CpG$Percent.Present)

# Add columns required by GrimAge clocks
Pheno$Female <- ifelse(Pheno$sex == "F", 1, 0)
Pheno$Age <- Pheno$age

# Define clocks to calculate
userClocks <- c("calcDNAmTL", "calcDamAge",
                "calcPhenoAge","calcHannum",
                 "calcHorvath2","calcHorvath1")

# Calculate all clocks stored in the vector
df_resultados <- calcUserClocks(
  userClocks,
  df_final,
  Pheno,
  imputation = FALSE
)

# Calculate both versions of GrimAge
GrimAgeV1 <- calcGrimAgeV1(df_final, Pheno, ID = "Sample_ID")
GrimAgeV2 <- calcGrimAgeV2(df_final, Pheno, ID = "Sample_ID")

# Join with the rest of results
df_resultados$GrimAgeV1 <- GrimAgeV1$GrimAgeV1
df_resultados$GrimAgeV2 <- GrimAgeV2$GrimAgeV2

# Get clock names (without 'calc' prefix)
# To identify columns in df_resultados
relojes <- c(gsub("calc", "", userClocks), "GrimAgeV1", "GrimAgeV2")


##### Generation of plots by separate tissue #####

df_resultados$Group <- as.character(df_resultados$Group)
df_resultados_blood <- filter(df_resultados, Group %in% c("control blood", "FA blood"))
df_resultados_bm <- filter(df_resultados, Group %in% c("control blood", "FA bone marrow"))

# Datasets dictionary
datasets <- list(
  Blood = list(df = df_resultados_blood),
  BM    = list(df = df_resultados_bm)
)

for (reloj in relojes) {
  cat("\n Procesando reloj:", reloj, "\n")
  
  plot_list <- list()
  
  for (name in names(datasets)) {
    # Get data
   dat <- datasets[[name]]$df
    
    try({
      # Scatter and linear regression plot
      p_reg <- ggplot(dat, aes(x = age, y = .data[[reloj]], color = Group)) +
        scale_color_manual(name = "Group", values = my_palette, labels = \(x) my_labels[x], drop = FALSE,
                           guide = "none")+
        geom_point(aes(fill = Group),
                   color = "black",
                   stroke = 0.2,
                   shape = 21,
                   size = 3) +
        scale_fill_manual(name = "Group", values = my_palette, labels = \(x) my_labels[x], drop = FALSE)+
        geom_smooth(aes(group = Group), method = "lm", se = FALSE, linewidth = 1.2) +
        # To add statistics:
        stat_cor(
          aes(label = paste(after_stat(r.label), after_stat(p.label), sep = "~`,`~")),
          label.x.npc = "center",      
          label.y.npc = "top",       
          size = 3,                  
          show.legend = FALSE        
        ) +
        scale_x_continuous(breaks = function(x) seq(0,100, by = 5)) +
        scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 5))+
        labs(title = paste(reloj, "(", name, ")"), x = "Age", y = "DNAm Age") +
        estilo_publicacion +
        theme(legend.position = "bottom",
              axis.title = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 14, color = "black"))
      
      # Residuals are already saved in the table
      nombre_res <- paste0("EAAres_", reloj) 
      # Statistic according to normality and homoscedasticity (T-test / Wilcoxon)
      formula_obj <- as.formula(paste(nombre_res, "~ Group"))
      
      # Normality test by group
      p_shapiro <- dat %>% 
        group_by(Group) %>% 
        shapiro_test(!!sym(nombre_res)) %>% 
        pull(p) %>% min(na.rm = TRUE)
      
      if (p_shapiro < 0.05) {
        pwc <- dat %>% wilcox_test(formula_obj) %>% add_significance()
        metodo_label <- "Wilcoxon"
      } else {
        # Parametric test
        # Check Variances (Levene)
        p_levene <- dat %>% levene_test(formula_obj) %>% pull(p)
        
        if (p_levene > 0.05) {
          # Equal variances: t test
        pwc <- dat %>% t_test(formula_obj, var.equal = TRUE) %>% add_significance()
        metodo_label <- "T-test"
        } else {
          # # Unequal variances: Welch's t-test
          pwc <- dat %>% t_test(formula_obj, var.equal = FALSE) %>% add_significance()
          metodo_label <- "Welch T-test"
        }}
      
      # Boxplot
      pwc <- pwc %>% add_xy_position(x = "Group")
      
      p_box <- ggplot(dat, aes(x = Group, y = .data[[nombre_res]], fill = Group)) +
        scale_fill_manual(name = "Group", values = my_palette, labels = \(x) my_labels[x], drop = FALSE,
                          guide = "none")+
        geom_boxplot(outlier.shape = NA, lwd = 0.2, median.linewidth= 0.5, width = 0.4, alpha = 0.2, show.legend = FALSE) +
        geom_quasirandom(aes(fill = Group), 
                         shape= 21,  color = "black",stroke = 0.1,
                         size = 2, width = 0.1, 
                         show.legend = FALSE) +
        scale_y_continuous(
         limits = c(NA, 50),                
          breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 10))+
        stat_summary(
          fun = median, 
          geom = "text", 
          aes(
            x = stage(Group, after_stat = as.numeric(x) + 0.3), 
            label = after_stat(round(y, 1))
          ), 
          vjust = 0.5,      
          hjust = 0,        
          color = "black", 
          fontface = "bold", 
          size = 3,
          inherit.aes = TRUE 
        ) +
        stat_pvalue_manual(pwc, label = "p",size = 3,  tip.length = 0.05) +
        scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1))+
        labs( y = "DNAm Age accel.", x = NULL, fill = "", subtitle = paste("Method:", metodo_label))+
        estilo_publicacion +
        theme(legend.position = "bottom",
              legend.key.size = unit(0.3, "cm"),
              axis.text.x = element_blank(), axis.ticks.x = element_blank())
      
      # Save to list for assembly
      plot_list[[paste0(name, "_reg")]] <- p_reg
      plot_list[[paste0(name, "_box")]] <- p_box
      
    }, silent = FALSE)
  }
  
  # Wrap blood and bone marrow in a single figure
    if (length(plot_list) == 4) {
      combined_plot <- plot_list[["Blood_reg"]] + plot_list[["Blood_box"]]+ 
        plot_list[["BM_reg"]] + plot_list[["BM_box"]] + 
        plot_layout(
          ncol = 2,
          widths = c(5, 2),
          guides = "keep" 
        ) 
    }
    
    # Save plots
    file_name <- paste0("Analysis_EAA_", reloj, ".png")
    ggsave(file.path(plot_path, file_name), 
           plot = combined_plot, width = 7, height = 10, dpi = 300)
    
    cat("Guardado:", file_name, "\n")
  }


# Save data objects

save(missing_CpG, df_final, Pheno, df_resultados, 
      df_resultados_blood, df_resultados_bm,  file = file.path(plot_path, "data_clocks.RData"), compress = TRUE)

