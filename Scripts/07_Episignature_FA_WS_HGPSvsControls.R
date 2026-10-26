##### Load R packages #####
library(ggplot2)
library(minfi)
library(limma)
library(tidyr)
library(dplyr)
library(ggpubr) 
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(readxl)
library(pheatmap)
library(tidyverse)
library(DMRcate) 
library(ggplotify)

##### Paths #####
path_base <- normalizePath(getwd())
input_path <- normalizePath(file.path(path_base, "Input_Betas"))
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))
signature_path <- normalizePath(file.path(path_base, "EpiSignature"))

##### Style configuration #####
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
##### Load results and prepare data frames #####

# Load results obtained from processing
beta_harmonized <- read_csv(file.path(input_path,"beta_completo.csv"))
Pheno <- read_csv(file.path(input_path,"all_targets_completo.csv"))

# Convert to data frame and set first column as row names
beta_harmonized <- as.data.frame(beta_harmonized) %>%
  column_to_rownames(var = "...1")
# Adjust Pheno data frame
Pheno$Sample_ID <- basename(Pheno$Basename) 
Pheno$Sample_ID <- ifelse(grepl("GSM", Pheno$Sample_ID), 
                          paste(Pheno$Sample_ID, Pheno$Condition), 
                          Pheno$Sample_ID)

Pheno$group <- paste(Pheno$Condition, Pheno$Tissue, sep = " ")

##### Obtaining the set of CpGs that compose the episignature #####
# Get annotation
ann_full <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Keep only probes with valid names present in both objects
common_probes <- intersect(rownames(beta_harmonized), rownames(ann_full))
common_probes <- common_probes[!is.na(common_probes)] 

beta_harmonized <- beta_harmonized[common_probes, ]

ann_filtered <- ann_full[common_probes, ]

# Create GRanges object specifying coordinate columns
gr_obj <- makeGRangesFromDataFrame(ann_filtered, 
                                   keep.extra.columns = TRUE,
                                   seqnames.field = "chr",
                                   start.field = "pos",
                                   end.field = "pos")

# Create GenomicRatioSet
grSet <- GenomicRatioSet(gr = gr_obj, 
                         Beta = as.matrix(beta_harmonized), 
                         annotation = c(array = "IlluminaHumanMethylationEPIC", 
                                        annotation = "ilm10b4.hg19"))

# Filter SNPs (MAF > 0.01) as mentioned in Pagliara et al.
grSet_clean <- dropLociWithSnps(grSet, snps = c("SBE", "CpG"), maf = 0.01)
beta_final_firma <- getBeta(grSet_clean)

# Load file of probes reported in Pagliara et al.
# Table S3. Probes defining the FANCA DNAm signature
# available in Supplemental information of the original article
probes_df <- read_excel(file.path(input_path,"probes_episignature_FA.xlsx"), 
                          sheet = "Table_probes_features")

# Extract probe IDs 
ids_episignature <- as.character(probes_df$ID)

# Check how many probes are present
sondas_presentes_clean <- intersect(rownames(beta_final_firma_clean), ids_episignature)

cat("Sondas de la epifirma encontradas:", length(sondas_presentes), "de", length(ids_episignature), "\n")
# We found 71 of 82

# Filter matrix to keep episignature probes
beta_firma_user <- beta_final_firma[sondas_presentes, ]

##### Multidimensional Scaling #####
# Use Euclidean distances as suggested by the article
dist_matrix <- dist(t(beta_firma_user))
mds_res <- cmdscale(dist_matrix)

# Create data frame for ggplot
df_plot_firma <- data.frame(
  MDS1 = mds_res[, 1],
  MDS2 = mds_res[, 2],
  Group = Pheno$group,
  Sex = Pheno$sex,
  Sample = Pheno$Sample_Name
)

# Generate plot
firma_mds <- ggplot(df_plot_firma, aes(x = MDS1, y = MDS2, fill = Group)) +
  geom_point(color = "black",
             stroke = 0.2,
             shape = 21,
             size = 4)+
  estilo_publicacion +
  scale_fill_manual(values = progeroid_palette, labels = progeroid_labels)+
  scale_y_continuous(expand = expansion(mult = 0.15)) + 
  scale_x_continuous(expand = expansion(mult = 0.15))+
 labs(
     x = "MDS 1",
    y = "MDS 2",
  )

ggsave((file.path(signature_path,"epifirma_mds.png")),
       plot = firma_mds,
       width = 7,
       height = 6,
       dpi = 300)

##### Heatmap #####
# Make sure Pheno table is in the same order as beta_firma_user columns
Pheno <- Pheno[match(colnames(beta_firma_user), Pheno$Sample_ID), ]

# Create annotation for heatmap columns
annotation_col <- data.frame(
  Group = factor(Pheno$group)
)
rownames(annotation_col) <- colnames(beta_firma_user)

ann_colors <- list(
  Group = progeroid_palette
)

# Plot Heatmap
firma_heatmap <- pheatmap(beta_firma_user,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2", 
         annotation_col = annotation_col,
         annotation_colors = list(
           Group = progeroid_palette),
         labels_col = Pheno$Sample_Name,
         show_colnames = FALSE, 
         fontsize_col = 6, 
         show_rownames = FALSE,
         main = paste("Fanconi Episignature Heatmap (", length(sondas_presentes), "probes)"),
         color = colorRampPalette(c("blue", "white", "yellow"))(100)) 


ggsave((file.path(signature_path,"epifirma_heatmap.jpg")),
plot = firma_heatmap,
width = 12,
height = 8,
dpi = 300)

