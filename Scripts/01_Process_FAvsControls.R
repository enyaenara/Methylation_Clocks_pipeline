##### Load R packages #####
library(minfi)
library(limma)
library(dplyr)
library(purrr)
library(tools)


##### Paths #####
path_base <- normalizePath(getwd())
data_path <- normalizePath(file.path(path_base, "Project_Data"))
idat_path <- normalizePath(file.path(path_base, "Raw_Data"))
rset_path <- normalizePath(file.path(path_base, "Rset_Data"))
input_path <- normalizePath(file.path(path_base, "Input_Betas/"))
normalized_path <- normalizePath(file.path(path_base, "Normalized_Data"))
results_path <- normalizePath(file.path(path_base, "Results"))
code_path <- normalizePath(file.path(path_base, "Scripts"))
plot_path <- normalizePath(file.path(path_base, "Plot_clocks"))

##### Functions to filter by control #####

getControlsPheno <- function(df, pheno_value) {
  
  
  df_control <- filter(df, Condition == pheno_value)
  
  
  return(df_control)
}

getControlsBetas <- function(df) {
  
  df_control <- dplyr::select(df, contains("control"))
  
  return(df_control)
}

##### Pre-run checks #####
# Search for physical files to ensure the path is real
for(prefix in c("GSE182991","GSE125367", #"GSE100825", discarded due to ages
                "GSE133774","GSE145233","GSE191276", 
                "Fanconi_anemia_blood", "Fanconi_anemia_bm",
                "Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm","ControlesEPICv2")) {
  
prefix_path <-  normalizePath(file.path(idat_path, prefix))
  
archivos_grn <- list.files(prefix_path, pattern = "_Grn.idat$", full.names = TRUE)
nombres_base <- gsub("_Grn.idat$", "", archivos_grn)
cat("Se encontraron todos los archivos de", prefix, "\n")
}

# Filter each folder to match the targets files
for(prefix in c("GSE182991","GSE125367",
                "GSE133774","GSE145233","GSE191276",
                "Fanconi_anemia_blood", "Fanconi_anemia_bm",
                "Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm","ControlesEPICv2")) {
  
  # Load targets and define paths
  targets_file <- file.path(data_path, paste0(prefix, "_targets.csv"))
  if(!file.exists(targets_file)) next 
  
  # We previously filter our csv files with age-range individuals
  targets <- read.csv(targets_file)
  prefix_path <- normalizePath(file.path(idat_path, prefix))
  older_path <- file.path(prefix_path, "older")
  
  # Build the filenames that should exist according to the table
  if (prefix %in% c("Fanconi_anemia_blood","Fanconi_anemia_bm",
                    "Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm","ControlesEPICv2")) {
    expected_basenames <- targets$Array_id
  } else {
    #For downloaded GEO files:
    expected_basenames <- paste(targets$Sample_Name, targets$Array_id, sep = "_")
  }
  
  # Create a list of complete files (Grn and Red)
  expected_files <- c(paste0(expected_basenames, "_Grn.idat"), 
                      paste0(expected_basenames, "_Red.idat"))
  
  # List actual files in the folder (only .idat)
  actual_files <- list.files(prefix_path, pattern = "\\.idat$")
  
  # Identify files that are not in the table 
  files_to_move <- setdiff(actual_files, expected_files)
  
  if (length(files_to_move) > 0) {
    if (!dir.exists(older_path)) dir.create(older_path)
    
    # Move leftover files
    file.rename(from = file.path(prefix_path, files_to_move),
                to = file.path(older_path, files_to_move))
    
    message(paste("En", prefix, "se movieron", length(files_to_move), "archivos a 'older'"))
  }
  
  # Update and save the Basename column in the target
  targets$Basename <- file.path(prefix_path, expected_basenames)
  write.csv(targets, file=file.path(prefix_path, "targets_match.csv"), row.names = FALSE)
  
  rm(targets)
  gc()
}

# Check that files exist before running
for(prefix in c("GSE182991","GSE125367",
                "GSE133774","GSE145233","GSE191276", 
                "Fanconi_anemia_blood","Fanconi_anemia_bm", 
                "Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm","ControlesEPICv2")) {

  prefix_path <-  normalizePath(file.path(idat_path, prefix))
  targets_final <- read.csv(file.path(prefix_path, "targets_match.csv"))
  test_file <- paste0(targets_final$Basename[1], "_Grn.idat")
if(file.exists(test_file)) {
  print("¡Éxito! R puede localizar los archivos IDAT.")
} else {
  stop("R aún no encuentra los archivos. Verifica que el nombre en 'Array_ID' coincida con el inicio del archivo IDAT.")
}
rm(targets_final)
gc()
}

##### Process EPIC v1 (Fanconi and GEO available datasets) #####
for(prefix in c("GSE182991","GSE125367",
                "GSE133774","GSE145233","GSE191276", 
                "Fanconi_anemia_blood",  "Fanconi_anemia_bm" )) {
  
  cat("Se esta procesando", prefix, "\n")
  prefix_path <- normalizePath(file.path(idat_path, prefix))
  targets_v1 <- read.csv(file.path(prefix_path, "targets_match.csv"))
  
  rgSet_v1_raw <- read.metharray.exp(targets = targets_v1,
                                 extended = TRUE,
                                 verbose = TRUE,
                                 force=TRUE)
  
 # Check probe quality 
  detP <- detectionP(rgSet_v1_raw)
  targets_v1$n_failed_probes <- colSums(detP > 0.01)
  targets_v1$mean_detection_p <- colMeans(detP)
 
  # Sample filtering (>5% failed)
  # to avoid introducing noise
  
  failsamp <- targets_v1$n_failed_probes > (0.05 * nrow(detP))
  
  # Create an rgSet without bad quality samples
  rgSet_v1 <- rgSet_v1_raw[, !failsamp]
  targets_v1 <- targets_v1[!failsamp, ]
  detP_clean <- detP[, !failsamp] # Also clean the p-values matrix

  # Save intermediate files 
   save(rgSet_v1_raw, rgSet_v1, detP_clean, failsamp,
       file = file.path(rset_path, prefix,
                        paste0(prefix, "_rset.RData")),
       compress = TRUE)
  

  # Normalization
  
  gSet_v1 <- preprocessFunnorm(rgSet_v1)
  
  
  # Probe filtering (Post-normalization)
  # Use 'keep' on the normalized object
  common_probes <- intersect(rownames(gSet_v1), rownames(detP_clean))
  keep <- rowSums(detP_clean[common_probes, ] < 0.01) == ncol(detP_clean)
  gSet_v1 <- gSet_v1[common_probes[keep], ]
  
  beta_v1 <- getBeta(gSet_v1)
  save(gSet_v1, targets_v1, beta_v1,
       file = file.path(normalized_path, prefix,
                        paste0(prefix, "_normalized.RData")),
       compress = TRUE)
  cat("Se proceso correctamente",prefix, "\n")
  rm(gSet_v1, targets_v1, beta_v1, detP,detP_clean,rgSet_v1,rgSet_v1_raw, keep, failsamp)
  gc()
}

# Merge all beta files into a single one to get the intersection of all EPICv1 datasets
for(prefix in c("Fanconi_anemia_blood", "Fanconi_anemia_bm","GSE182991","GSE125367",
                "GSE133774","GSE145233","GSE191276")) {
  load(file=file.path(normalized_path, prefix, paste(prefix, "normalized.RData", sep="_")))
  if (prefix %in% c("Fanconi_anemia_blood", "Fanconi_anemia_bm")) {
  # Save targets and betas table of Fanconi anemia
    targets<- targets_v1
    new_name <- paste(prefix, "targets", sep = "_" )
    assign(new_name, targets)
    
    betas <- beta_v1
    new_name <- paste(prefix,"betas", sep = "_" )
    assign(new_name, betas)
    
    } else {
  # Get all control targets and name with prefix
  targets<- getControlsPheno(targets_v1, "control")
  new_name <- paste(prefix, "targets", sep = "_" )
  assign(new_name, targets)
  
  # Get all control betas and name with prefix
  pheno <- targets_v1$Condition
  colnames(beta_v1) <- paste(colnames(beta_v1), pheno)
  betas <- getControlsBetas(as.data.frame(beta_v1))
  new_name <- paste(prefix,"betas", sep = "_" )
  assign(new_name, betas)
   }
  rm(beta_v1,betas, targets, gSet_v1,targets_v1)
  gc()
   
}

# Bind target rows
x <- list(Fanconi_anemia_blood_targets, Fanconi_anemia_bm_targets,
          GSE125367_targets,
          GSE133774_targets,GSE145233_targets,
          GSE182991_targets,GSE191276_targets)

all_targets_v1 <- do.call(rbind, x)
all_targets_v1$age <- round(as.numeric(all_targets_v1$age), 1)

# Remove files no longer needed
rm(Fanconi_anemia_blood_targets, Fanconi_anemia_bm_targets,
   GSE125367_targets,
   GSE133774_targets,GSE145233_targets,
   GSE182991_targets,GSE191276_targets)

# Create CpG column to join into a single data frame
x <- list(Fanconi_anemia_blood_betas, Fanconi_anemia_bm_betas,
          GSE125367_betas,
          GSE133774_betas,GSE145233_betas,
          GSE182991_betas,GSE191276_betas)

common_cpgs <- Reduce(intersect, lapply(x, rownames))

all_betas_v1 <- do.call(cbind, lapply(x, function(df) {
  df[common_cpgs, , drop = FALSE]
}))

rm(Fanconi_anemia_blood_betas, Fanconi_anemia_bm_betas,
   GSE125367_betas,
   GSE133774_betas,GSE145233_betas,
   GSE182991_betas,GSE191276_betas)
gc()

##### Process EPIC v2 (Fanconi and Mexican controls) #####
for(prefix in c("Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm", "ControlesEPICv2")) {
  cat("Se esta procesando", prefix, "\n")
  prefix_path <- normalizePath(file.path(idat_path, prefix))
targets_v2 <- read.csv(file.path(prefix_path, "targets_match.csv"))
rgSet_v2_raw <- read.metharray.exp(targets = targets_v2,
                               extended = TRUE,
                               verbose = TRUE)


# Assign annotation so minfi recognizes v2
annotation(rgSet_v2_raw) <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")

# Check probe quality
detP <- detectionP(rgSet_v2_raw)
targets_v2$n_failed_probes <- colSums(detP > 0.01)
targets_v2$mean_detection_p <- colMeans(detP)

# Sample filtering (>5% failed)

failsamp <- targets_v2$n_failed_probes > (0.05 * nrow(detP))

rgSet_v2 <- rgSet_v2_raw[, !failsamp]
targets_v2 <- targets_v2[!failsamp, ]
detP_clean <- detP[, !failsamp] 

# Save intermediate files
save(rgSet_v2_raw, rgSet_v2, detP, detP_clean, failsamp,
     file = file.path(rset_path, prefix,
                      paste0(prefix, "_rset.RData")),
     compress = TRUE)


# Normalize v2
gSet_v2 <- preprocessFunnorm(rgSet_v2)

# Probe filtering (Post-normalization)

common_probes <- intersect(rownames(gSet_v2), rownames(detP_clean))
keep <- rowSums(detP_clean[common_probes, ] < 0.01) == ncol(detP_clean)
gSet_v2 <- gSet_v2[common_probes[keep], ]


beta_v2 <- getBeta(gSet_v2)

# Remove suffixes (_BCxx) so names match v1
rownames(beta_v2) <- gsub("_.*$", "", rownames(beta_v2))

# Harmonization (v2 contains multiple probes for the same CpG)
# This averages replicates and leaves a single row for each cgXXXXXXXX
beta_v2_clean <- avearrays(beta_v2)


save(gSet_v2, targets_v2, beta_v2,beta_v2_clean,
     file = file.path(normalized_path, prefix,
                      paste0(prefix, "_normalized.RData")),
     compress = TRUE)

cat("Se proceso correctamente",prefix, "\n")
rm(gSet_v2, targets_v2, beta_v2, beta_v2_clean, detP,detP_clean,rgSet_v2,rgSet_v2_raw, keep, failsamp)
gc()
}

# Merge all beta files into a single one to get the intersection of all EPICv2 datasets
for(prefix in c("Fanconi_anemia_V2_blood", "Fanconi_anemia_V2_bm", "ControlesEPICv2")) {
  load(file=file.path(normalized_path, prefix, paste(prefix, "normalized.RData", sep="_")))
  
  # Save targets and betas matrices
  targets<- targets_v2
  new_name <- paste(prefix, "targets", sep = "_" )
  assign(new_name, targets)
  
  betas <- beta_v2_clean
  new_name <- paste(prefix,"betas", sep = "_" )
  assign(new_name, betas)
  
  
  rm(beta_v2, beta_v2_clean, betas, targets, gSet_v2,targets_v2)
  gc()
  
}

# Bind target rows
x <- list(Fanconi_anemia_V2_blood_targets, Fanconi_anemia_V2_bm_targets,
          ControlesEPICv2_targets)

all_targets_v2 <- do.call(rbind, x)
all_targets_v2$age <- round(as.numeric(all_targets_v2$age), 1)

# Remove files no longer needed
rm(Fanconi_anemia_V2_blood_targets, Fanconi_anemia_V2_bm_targets,
   ControlesEPICv2_targets)

# Create CpG column to join into a single data frame
x <- list(Fanconi_anemia_V2_blood_betas, Fanconi_anemia_V2_bm_betas,
          ControlesEPICv2_betas)

common_cpgs <- Reduce(intersect, lapply(x, rownames))

all_betas_v2_clean <- do.call(cbind, lapply(x, function(df) {
  df[common_cpgs, , drop = FALSE]
}))

rm(Fanconi_anemia_V2_blood_betas, Fanconi_anemia_V2_bm_betas,
   ControlesEPICv2_betas)
gc()


##### Final merge: Intersection of common probes ##### 
# Look for sites present in both (~700,00 sites)
common_cpgs <- intersect(rownames(all_betas_v1), rownames(all_betas_v2_clean))
print(paste("Sondas en común encontradas:", length(common_cpgs)))

# Merge matrices
beta_combined <- cbind(all_betas_v1[common_cpgs, ], all_betas_v2_clean[common_cpgs, ])

# Save betas for clock analysis
write.csv(beta_combined, file.path(input_path, "ALL_FA_vs_ALL/beta_harmonized.csv"), row.names = TRUE)

# Save data
targets_final <- rbind(all_targets_v1, all_targets_v2)
write.csv(targets_final,file = file.path(input_path, "ALL_FA_vs_ALL/all_targets.csv"), row.names=F)



