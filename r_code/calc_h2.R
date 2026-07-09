install.packages("tidyverse")
install.packages("gifti")
install.packages("solarius")
install.packages("ggpubr")
install.packages("Cairo")
install.packages("devtools")
#install_github("ugcd/solarius")
install.packages("R.matlab")
install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(doParallel)
library(devtools)
library(tidyverse)
library(gifti)
library(R.matlab)
library(solarius)
library(ggpubr)
library(Cairo)

# --- USER CONFIGURATION: edit this to point at your data directory ---
base_dir <- "/path/to/isc_heritability_r"
data_dir <- "/path/to/data"

num_days  <- 2    # scan days/repeats (Day 1 / Day 2)
num_edges <- 153  # unique network-pair edges (17 Kong/Yeo networks -> 17*18/2 = 153)

movie_rest_subj_indices <- c(1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,175,176,177,178)
perm_start <- 1
# --------------------------------------------------------------------

num_cores <- detectCores()  # Detect the number of cores
cl <- makeCluster(num_cores-1)  # Create a cluster
registerDoParallel(cl)  # Register the cluster



calc_net_heritability <- function(df, trait_list, covars, dir_string) {
  covar_formula <- paste0(covars, collapse = ' + ')
  h2_df <- NULL
  expected_cols <- NULL  # Variable to store the expected number of columns
  
  for (trait in trait_list) {
    write(trait, '')
    
    tryCatch({
      df[[trait]] <- scale(df[[trait]])[, 1]
      formula <- as.formula(paste0(trait, ' ~ ', covar_formula))
      rhog <- invisible(solarPolygenic(formula, df, dir = dir_string))
      
      out_row <- rhog$vcf[1,]
      print(out_row)
      out_row$trait <- trait
      
      # If h2_df is NULL, initialize it with the first out_row and store the expected number of columns
      if (is.null(h2_df)) {
        h2_df <- out_row
        expected_cols <- ncol(h2_df)
      } else {
        # Check if out_row has the same number of columns as expected
        if (ncol(out_row) == expected_cols) {
          h2_df <- rbind(h2_df, out_row)
        } else {
          # If not, create a row of NaNs with the correct number of columns
          nan_row <- as.data.frame(matrix(NaN, nrow = 1, ncol = expected_cols))
          colnames(nan_row) <- colnames(h2_df)
          nan_row$trait <- trait
          h2_df <- rbind(h2_df, nan_row)
        }
      }
    }, error = function(e) {
      # In case of an error, also return a row of NaNs
      warning(paste("Error in processing trait:", trait, "; Error:", e$message))
      if (!is.null(expected_cols)) {
        nan_row <- as.data.frame(matrix(NaN, nrow = 1, ncol = expected_cols))
        colnames(nan_row) <- colnames(h2_df)
        nan_row$trait <- trait
        h2_df <- rbind(h2_df, nan_row)
      }
    })
  }
  return(h2_df)
}



# Connectivity/Movie FC Magnitude
all_results <- list()
align_type <- "connectivity"
parc_list <- as.character(seq(300, 1000, by = 100))
task_type <- "movie"
do_jackknife <-0
for (parc in parc_list) {
  for (scan_id in 1:num_days) {
    data_df <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
    data_df <- as.data.frame(t(data_df))
    # Assuming 'data_df' is your dataframe
    zscored_data_df <- as.data.frame(scale(data_df))
    
    colnames(zscored_data_df) <- paste0("net", 1:ncol(data_df))
    colnames(data_df) <- paste0("net", 1:ncol(data_df))
    
    covari_df <- read.csv(paste0(base_dir, "/covariates_connectivity_scan_", scan_id, ".csv"), header = FALSE)
    colnames(covari_df) <- c("age", "sex", "motion")
    
    kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_rest.csv"))
    combined_df <- cbind(kinship_df, data_df, covari_df)
    
    # Ensure that FAMID is treated as a character
    combined_df$FAMID <- as.character(combined_df$FAMID)
    
    # Create a mapping from the unique family IDs to a sequence of integers
    unique_famids <- unique(combined_df$FAMID)
    famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
    
    # Apply the mapping to the FAMID column
    combined_df$FAMID <- famid_mapping[combined_df$FAMID]
    
    # Exclude NaN values
    # Count the number of occurrences for each value
    counts <- table(combined_df$MZTWIN)
    # Identify unique values (occurring only once)
    unique_values <- names(counts)[counts == 1]
    # Retrieve corresponding subject IDs
    combined_df$MZTWIN[combined_df$MZTWIN %in% unique_values] <- NaN
    dir_string = paste0(base_dir, "/fc")
    
    # Run heritability calculation with all families
    results_full <- calc_net_heritability(combined_df, colnames(data_df), colnames(covari_df),dir_string)
    all_results[[paste0("full_", scan_id)]] <- results_full
    write.csv(results_full, paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), row.names = FALSE)
    
    num_fams = length(unique_famids)
    if (do_jackknife == 1) {
      # Parallel Jackknife procedure
      foreach(fam_id = seq(num_fams), .combine = 'c', .packages = c('solarius')) %dopar% {
        tryCatch({
          save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
          if (!file.exists(save_path)) {
            df_jackknife <- combined_df[combined_df$FAMID != fam_id, ]
            dir_string = paste0(base_dir, "/fc/",fam_id)
            results_jackknife <- calc_net_heritability(df_jackknife, colnames(data_df), colnames(covari_df),dir_string)
            
            # Save jackknife results
            save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
            write.csv(results_jackknife, save_path, row.names = FALSE)
          } else {
          }
        }, error = function(e) {
          NA
        })
      }
    }
  }
}

# Piecewise/Movie FC Magnitude
all_results <- list()
align_type <- "piecewise"
parc_list <- as.character(seq(100, 1000, by = 100))
task_type <- "movie"
do_jackknife <- 0
for (parc in parc_list) {
  for (scan_id in 1:num_days) {
    data_df <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
    data_df <- as.data.frame(t(data_df))
    # Assuming 'data_df' is your dataframe
    zscored_data_df <- as.data.frame(scale(data_df))
    
    colnames(zscored_data_df) <- paste0("net", 1:ncol(data_df))
    colnames(data_df) <- paste0("net", 1:ncol(data_df))
    
    covari_df <- read.csv(paste0(base_dir, "/covariates_",task_type,"_scan_", scan_id, ".csv"), header = FALSE)
    colnames(covari_df) <- c("age", "sex", "motion")
    
    kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_",task_type,".csv"))
    combined_df <- cbind(kinship_df, data_df, covari_df)
    
    # Ensure that FAMID is treated as a character
    combined_df$FAMID <- as.character(combined_df$FAMID)
    
    # Create a mapping from the unique family IDs to a sequence of integers
    unique_famids <- unique(combined_df$FAMID)
    famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
    
    # Apply the mapping to the FAMID column
    combined_df$FAMID <- famid_mapping[combined_df$FAMID]
    
    # Exclude NaN values
    # Count the number of occurrences for each value
    counts <- table(combined_df$MZTWIN)
    # Identify unique values (occurring only once)
    unique_values <- names(counts)[counts == 1]
    # Retrieve corresponding subject IDs
    combined_df$MZTWIN[combined_df$MZTWIN %in% unique_values] <- NaN
    dir_string = paste0(base_dir, "/fc")
    
    # Run heritability calculation with all families
    results_full <- calc_net_heritability(combined_df, colnames(data_df), colnames(covari_df),dir_string)
    all_results[[paste0("full_", scan_id)]] <- results_full
    write.csv(results_full, paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), row.names = FALSE)
    if (do_jackknife == 1) {
      
      num_fams = length(unique_famids)
      # Parallel Jackknife procedure
      
      foreach(fam_id = seq(num_fams), .combine = 'c', .packages = c('solarius')) %dopar% {
        tryCatch({
          save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
          if (!file.exists(save_path)) {
            df_jackknife <- combined_df[combined_df$FAMID != fam_id, ]
            dir_string = paste0(base_dir, "/fc/",fam_id)
            results_jackknife <- calc_net_heritability(df_jackknife, colnames(data_df), colnames(covari_df),dir_string)
            
            # Save jackknife results
            save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
            write.csv(results_jackknife, save_path, row.names = FALSE)
          } else {
          }
        }, error = function(e) {
          NA
        })
      }
    }
  }
}

# Anatomical/Movie FC Magnitude
all_results <- list()
align_type <- "anatomical"
parc_list <- as.character(seq(400, 400, by = 100))
task_type <- "movie"
do_jackknife <- 0
for (parc in parc_list) {
  for (scan_id in 1:num_days) {
    data_df <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
    data_df <- as.data.frame(t(data_df))
    # Assuming 'data_df' is your dataframe
    zscored_data_df <- as.data.frame(scale(data_df))
    
    colnames(zscored_data_df) <- paste0("net", 1:ncol(data_df))
    colnames(data_df) <- paste0("net", 1:ncol(data_df))
    
    covari_df <- read.csv(paste0(base_dir, "/covariates_",task_type,"_scan_", scan_id, ".csv"), header = FALSE)
    colnames(covari_df) <- c("age", "sex", "motion")
    
    kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_",task_type,".csv"))
    combined_df <- cbind(kinship_df, data_df, covari_df)
    
    # Ensure that FAMID is treated as a character
    combined_df$FAMID <- as.character(combined_df$FAMID)
    
    # Create a mapping from the unique family IDs to a sequence of integers
    unique_famids <- unique(combined_df$FAMID)
    famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
    
    # Apply the mapping to the FAMID column
    combined_df$FAMID <- famid_mapping[combined_df$FAMID]
    
    # Exclude NaN values
    # Count the number of occurrences for each value
    counts <- table(combined_df$MZTWIN)
    # Identify unique values (occurring only once)
    unique_values <- names(counts)[counts == 1]
    # Retrieve corresponding subject IDs
    combined_df$MZTWIN[combined_df$MZTWIN %in% unique_values] <- NaN
    dir_string = paste0(base_dir, "/fc")
    
    # Run heritability calculation with all families
    results_full <- calc_net_heritability(combined_df, colnames(data_df), colnames(covari_df),dir_string)
    all_results[[paste0("full_", scan_id)]] <- results_full
    write.csv(results_full, paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), row.names = FALSE)
    if (do_jackknife == 1) {
      
      num_fams = length(unique_famids)
      # Parallel Jackknife procedure
      
      foreach(fam_id = seq(num_fams), .combine = 'c', .packages = c('solarius')) %dopar% {
        df_jackknife <- combined_df[combined_df$FAMID != fam_id, ]
        dir_string = paste0(base_dir, "/fc/",fam_id)
        results_jackknife <- calc_net_heritability(df_jackknife, colnames(data_df), colnames(covari_df),dir_string)
        
        # Save jackknife results
        save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
        write.csv(results_jackknife, save_path, row.names = FALSE)
      }
    }
  }
}

# Anatomical/Rest FC magnitude
all_results <- list()
align_type <- "anatomical"
parc_list <- as.character(seq(400, 400, by = 100))
task_type <- "rest"
do_jackknife <- 0
for (parc in parc_list) {
  for (scan_id in 1:num_days) {
    data_df <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
    data_df <- as.data.frame(t(data_df))
    # Assuming 'data_df' is your dataframe
    zscored_data_df <- as.data.frame(scale(data_df))
    
    colnames(zscored_data_df) <- paste0("net", 1:ncol(data_df))
    colnames(data_df) <- paste0("net", 1:ncol(data_df))
    
    covari_df <- read.csv(paste0(base_dir, "/covariates_",task_type,"_scan_", scan_id, ".csv"), header = FALSE)
    colnames(covari_df) <- c("age", "sex", "motion")
    
    kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_",task_type,".csv"))
    combined_df <- cbind(kinship_df, data_df, covari_df)
    
    # Ensure that FAMID is treated as a character
    combined_df$FAMID <- as.character(combined_df$FAMID)
    
    # Create a mapping from the unique family IDs to a sequence of integers
    unique_famids <- unique(combined_df$FAMID)
    famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
    
    # Apply the mapping to the FAMID column
    combined_df$FAMID <- famid_mapping[combined_df$FAMID]
    
    # Exclude NaN values
    # Count the number of occurrences for each value
    counts <- table(combined_df$MZTWIN)
    # Identify unique values (occurring only once)
    unique_values <- names(counts)[counts == 1]
    # Retrieve corresponding subject IDs
    combined_df$MZTWIN[combined_df$MZTWIN %in% unique_values] <- NaN
    
    # Run heritability calculation with all families
    dir_string = paste0(base_dir, "/fc")
    
    results_full <- calc_net_heritability(combined_df, colnames(data_df), colnames(covari_df),dir_string)
    all_results[[paste0("full_", scan_id)]] <- results_full
    write.csv(results_full, paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), row.names = FALSE)
    if (do_jackknife == 1) {
      num_fams = length(unique_famids)
      # Parallel Jackknife procedure
      foreach(fam_id = seq(num_fams), .combine = 'c', .packages = c('solarius')) %dopar% {
        df_jackknife <- combined_df[combined_df$FAMID != fam_id, ]
        dir_string = paste0(base_dir, "/fc/",fam_id)
        results_jackknife <- calc_net_heritability(df_jackknife, colnames(data_df), colnames(covari_df),dir_string)
        
        # Save jackknife results
        save_path <- paste0(data_dir, "/isc_heritability/data/solar/fc/fc_herit_net_",task_type,"_",align_type,"_parc_",parc,"_scan_",scan_id, "_jackknife_", fam_id, ".csv")
        write.csv(results_jackknife, save_path, row.names = FALSE)
      }
    }
  }
}

############
# Basic setup

# Load the required libraries
library(foreach)
library(doParallel)

# The permutation section below is far more compute-intensive than the code
# above, so it gets its own (typically larger) cluster. Stop the earlier
# cluster first so both aren't left running simultaneously.
stopCluster(cl)

# Number of cores to use for the permutation section (edit for your machine)
num_cores_perm <- 16

# Register parallel backend
cl <- makeCluster(num_cores_perm)
registerDoParallel(cl)
clusterExport(cl, c("calc_net_heritability"))

# Number of permutations
num_perm <- 1000
############
# Generate null differences for comparisons of interest


# Movie vs. Rest
all_results <- list()
align_type <- "anatomical"
parc_list <- '400'
task_type1 <- "rest"
task_type2 <- "movie"
do_jackknife <-0

# Parallelized loop
x<-foreach(perm_id = 1:num_perm, .packages = c("data.table","solarius")) %dopar% {
  print(perm_id)
  for (parc in parc_list) {
    for (scan_id in 1:num_days) {
      data_df1 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type1,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
      data_df1 <- as.data.frame(t(data_df1))
      data_df2 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type2,"_",align_type,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
      data_df2 <- as.data.frame(t(data_df2))
      data_df2 <- data_df2[unlist(movie_rest_subj_indices), ]
      data_df1 <- as.data.frame(scale(data_df1))
      data_df2 <- as.data.frame(scale(data_df2))
      colnames(data_df1) <- paste0("net", 1:ncol(data_df1))
      colnames(data_df2) <- paste0("net", 1:ncol(data_df2))
      covari_df1 <- read.csv(paste0(base_dir, "/covariates_",task_type1,"_scan_", scan_id, ".csv"), header = FALSE)
      colnames(covari_df1) <- c("age", "sex", "motion")
      covari_df2 <- read.csv(paste0(base_dir, "/covariates_",task_type2,"_scan_", scan_id, ".csv"), header = FALSE)
      colnames(covari_df2) <- c("age", "sex", "motion")
      covari_df2 <- covari_df2[unlist(movie_rest_subj_indices), ]
      
      num_subj_swap <- nrow(data_df1)
      set.seed(perm_id)
      
      # Generate random array of 1s and 2s
      swap_indices <- sample(1:2, num_subj_swap, replace = TRUE)
      
      # Initialize swapped data frames
      swapped_data_df1 <- data_df1
      swapped_data_df2 <- data_df2
      swapped_covari_df1 <- covari_df1
      swapped_covari_df2 <- covari_df2
      
      # Swap rows between data frames based on swap_indices
      for (i in 1:num_subj_swap) {
        if (swap_indices[i] == 1) {
          temp_row <- swapped_data_df1[i, ]
          swapped_data_df1[i, ] <- swapped_data_df2[i, ]
          swapped_data_df2[i, ] <- temp_row
          
          temp_row_covari <- swapped_covari_df1[i, ]
          swapped_covari_df1[i, ] <- swapped_covari_df2[i, ]
          swapped_covari_df2[i, ] <- temp_row_covari
        }
      }
      
      kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_",task_type1,".csv"))
      #kinship_df <- kinship_df[unlist(movie_rest_subj_indices), ]
      
      combined_df1 <- cbind(kinship_df, swapped_data_df1, covari_df1)
      combined_df2 <- cbind(kinship_df, swapped_data_df2, covari_df2)
      
      # Ensure that FAMID is treated as a character
      combined_df1$FAMID <- as.character(combined_df1$FAMID)
      combined_df2$FAMID <- as.character(combined_df2$FAMID)
      
      # Create a mapping from the unique family IDs to a sequence of integers
      unique_famids <- unique(combined_df1$FAMID)
      famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
      
      # Apply the mapping to the FAMID column
      combined_df1$FAMID <- famid_mapping[combined_df1$FAMID]
      combined_df2$FAMID <- famid_mapping[combined_df2$FAMID]
      
      # Exclude NaN values
      # Count the number of occurrences for each value
      counts <- table(combined_df1$MZTWIN)
      # Identify unique values (occurring only once)
      unique_values <- names(counts)[counts == 1]
      # Retrieve corresponding subject IDs
      combined_df1$MZTWIN[combined_df1$MZTWIN %in% unique_values] <- NaN
      combined_df2$MZTWIN[combined_df2$MZTWIN %in% unique_values] <- NaN
      
      dir_string1 <- paste0(base_dir, "/fc_perm/perm_", perm_id, "/")
      dir_string2 <- paste0(base_dir, "/fc_perm2/perm_", perm_id, "/")
      
      # Create directories if they don't exist
      if (!file.exists(dir_string1)) {
        dir.create(dir_string1, recursive = TRUE)
      }
      if (!file.exists(dir_string2)) {
        dir.create(dir_string2, recursive = TRUE)
      }
      
      # Your existing code here...
      
      results_full_1 <- calc_net_heritability(combined_df1, colnames(data_df1), colnames(covari_df1), dir_string1)
      results_full_2 <- calc_net_heritability(combined_df2, colnames(data_df2), colnames(covari_df2), dir_string2)
      results_diff1 <- results_full_2$Var - results_full_1$Var
      
      saveRDS(results_diff1, file = paste0(base_dir, "/fc_perm/perm_", perm_id, "/results_diff_scan_", scan_id, ".rds"))
    }
  }
}


# Piecewise vs. Anatomical
all_results <- list()
align_type1 <- "anatomical"
align_type2 <- "piecewise"
parc_list <- '400'
task_type <- "movie"
do_jackknife <-0

x<-foreach(perm_id = 1:num_perm, .packages = c("data.table","solarius")) %dopar% {
  print(perm_id)
  #for (perm_id in 1:num_perm) {
  for (parc in parc_list) {
    for (scan_id in 1:num_days) {
      data_df1 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type1,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
      data_df1 <- as.data.frame(t(data_df1))
      data_df2 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type2,"_parc_","100","_scan_",scan_id, ".csv"), header = FALSE)
      data_df2 <- as.data.frame(t(data_df2))
      data_df1 <- as.data.frame(scale(data_df1))
      data_df2 <- as.data.frame(scale(data_df2))
      colnames(data_df1) <- paste0("net", 1:ncol(data_df1))
      colnames(data_df2) <- paste0("net", 1:ncol(data_df2))
      covari_df1 <- read.csv(paste0(base_dir, "/covariates_",task_type,"_scan_", scan_id, ".csv"), header = FALSE)
      colnames(covari_df1) <- c("age", "sex", "motion")
      covari_df2 = covari_df1
      
      num_subj_swap <- nrow(data_df1)
      set.seed(perm_id)
      
      # Generate random array of 1s and 2s
      swap_indices <- sample(1:2, num_subj_swap, replace = TRUE)
      
      # Initialize swapped data frames
      swapped_data_df1 <- data_df1
      swapped_data_df2 <- data_df2
      swapped_covari_df1 <- covari_df1
      swapped_covari_df2 <- covari_df2
      
      # Swap rows between data frames based on swap_indices
      for (i in 1:num_subj_swap) {
        if (swap_indices[i] == 1) {
          temp_row <- swapped_data_df1[i, ]
          swapped_data_df1[i, ] <- swapped_data_df2[i, ]
          swapped_data_df2[i, ] <- temp_row
          
          temp_row_covari <- swapped_covari_df1[i, ]
          swapped_covari_df1[i, ] <- swapped_covari_df2[i, ]
          swapped_covari_df2[i, ] <- temp_row_covari
        }
      }
      
      kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_",task_type,".csv"))
      #kinship_df <- kinship_df[unlist(movie_rest_subj_indices), ]
      
      combined_df1 <- cbind(kinship_df, swapped_data_df1, covari_df1)
      combined_df2 <- cbind(kinship_df, swapped_data_df2, covari_df2)
      
      # Ensure that FAMID is treated as a character
      combined_df1$FAMID <- as.character(combined_df1$FAMID)
      combined_df2$FAMID <- as.character(combined_df2$FAMID)
      
      # Create a mapping from the unique family IDs to a sequence of integers
      unique_famids <- unique(combined_df1$FAMID)
      famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
      
      # Apply the mapping to the FAMID column
      combined_df1$FAMID <- famid_mapping[combined_df1$FAMID]
      combined_df2$FAMID <- famid_mapping[combined_df2$FAMID]
      
      # Exclude NaN values
      # Count the number of occurrences for each value
      counts <- table(combined_df1$MZTWIN)
      # Identify unique values (occurring only once)
      unique_values <- names(counts)[counts == 1]
      # Retrieve corresponding subject IDs
      combined_df1$MZTWIN[combined_df1$MZTWIN %in% unique_values] <- NaN
      combined_df2$MZTWIN[combined_df2$MZTWIN %in% unique_values] <- NaN
      
      dir_string1 <- paste0(base_dir, "/fc_perm_piecewise1/perm_", perm_id, "/")
      dir_string2 <- paste0(base_dir, "/fc_perm_piecewise2/perm_", perm_id, "/")
      
      # Create directories if they don't exist
      if (!file.exists(dir_string1)) {
        dir.create(dir_string1, recursive = TRUE)
      }
      if (!file.exists(dir_string2)) {
        dir.create(dir_string2, recursive = TRUE)
      }
      
      # Your existing code here...
      
      results_full_1 <- calc_net_heritability(combined_df1, colnames(data_df1), colnames(covari_df1), dir_string1)
      results_full_2 <- calc_net_heritability(combined_df2, colnames(data_df2), colnames(covari_df2), dir_string2)
      results_diff1 <- results_full_2$Var - results_full_1$Var
      
      saveRDS(results_diff1, file = paste0(base_dir, "/fc_perm_piecewise1/perm_", perm_id, "/results_diff_scan_", scan_id, ".rds"))
      
    }
  }
}

# Connectivity vs. Anatomical
all_results <- list()
align_type1 <- "anatomical"
align_type2 <- "connectivity"
parc_list <- '400'
task_type <- "movie"
do_jackknife <-0

x<-foreach(perm_id = perm_start:num_perm, .packages = c("data.table","solarius")) %dopar% {
  print(perm_id)
  #for (perm_id in 1:num_perm) {
  for (parc in parc_list) {
    for (scan_id in 1:num_days) {
      data_df1 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type1,"_parc_",parc,"_scan_",scan_id, ".csv"), header = FALSE)
      data_df1 <- as.data.frame(t(data_df1))
      data_df2 <- read.csv(paste0(base_dir, "/fc/fc_net_",task_type,"_",align_type2,"_parc_","100","_scan_",scan_id, ".csv"), header = FALSE)
      data_df2 <- as.data.frame(t(data_df2))
      data_df1 <- as.data.frame(scale(data_df1))
      data_df2 <- as.data.frame(scale(data_df2))
      data_df1 <- data_df1[unlist(movie_rest_subj_indices), ]
      
      colnames(data_df1) <- paste0("net", 1:ncol(data_df1))
      colnames(data_df2) <- paste0("net", 1:ncol(data_df2))
      covari_df1 <- read.csv(paste0(base_dir, "/covariates_connectivity_scan_", scan_id, ".csv"), header = FALSE)
      colnames(covari_df1) <- c("age", "sex", "motion")
      covari_df2 = covari_df1
      
      num_subj_swap <- nrow(data_df1)
      set.seed(perm_id)
      
      # Generate random array of 1s and 2s
      swap_indices <- sample(1:2, num_subj_swap, replace = TRUE)
      
      # Initialize swapped data frames
      swapped_data_df1 <- data_df1
      swapped_data_df2 <- data_df2
      swapped_covari_df1 <- covari_df1
      swapped_covari_df2 <- covari_df2
      
      # Swap rows between data frames based on swap_indices
      for (i in 1:num_subj_swap) {
        if (swap_indices[i] == 1) {
          temp_row <- swapped_data_df1[i, ]
          swapped_data_df1[i, ] <- swapped_data_df2[i, ]
          swapped_data_df2[i, ] <- temp_row
          
          temp_row_covari <- swapped_covari_df1[i, ]
          swapped_covari_df1[i, ] <- swapped_covari_df2[i, ]
          swapped_covari_df2[i, ] <- temp_row_covari
        }
      }
      
      kinship_df <- read.csv(paste0(data_dir, "/isc_heritability/pheno_table_rest.csv"))
      #kinship_df <- kinship_df[unlist(movie_rest_subj_indices), ]
      
      combined_df1 <- cbind(kinship_df, swapped_data_df1, covari_df1)
      combined_df2 <- cbind(kinship_df, swapped_data_df2, covari_df2)
      
      # Ensure that FAMID is treated as a character
      combined_df1$FAMID <- as.character(combined_df1$FAMID)
      combined_df2$FAMID <- as.character(combined_df2$FAMID)
      
      # Create a mapping from the unique family IDs to a sequence of integers
      unique_famids <- unique(combined_df1$FAMID)
      famid_mapping <- setNames(seq_along(unique_famids), unique_famids)
      
      # Apply the mapping to the FAMID column
      combined_df1$FAMID <- famid_mapping[combined_df1$FAMID]
      combined_df2$FAMID <- famid_mapping[combined_df2$FAMID]
      
      # Exclude NaN values
      # Count the number of occurrences for each value
      counts <- table(combined_df1$MZTWIN)
      # Identify unique values (occurring only once)
      unique_values <- names(counts)[counts == 1]
      # Retrieve corresponding subject IDs
      combined_df1$MZTWIN[combined_df1$MZTWIN %in% unique_values] <- NaN
      combined_df2$MZTWIN[combined_df2$MZTWIN %in% unique_values] <- NaN
      
      dir_string1 <- paste0(base_dir, "/fc_perm_connectivity1/perm_", perm_id, "/")
      dir_string2 <- paste0(base_dir, "/fc_perm_connectivity2/perm_", perm_id, "/")
      
      # Create directories if they don't exist
      if (!file.exists(dir_string1)) {
        dir.create(dir_string1, recursive = TRUE)
      }
      if (!file.exists(dir_string2)) {
        dir.create(dir_string2, recursive = TRUE)
      }
      
      # Your existing code here...
      
      results_full_1 <- calc_net_heritability(combined_df1, colnames(data_df1), colnames(covari_df1), dir_string1)
      results_full_2 <- calc_net_heritability(combined_df2, colnames(data_df2), colnames(covari_df2), dir_string2)
      results_diff1 <- results_full_2$Var - results_full_1$Var
      
      saveRDS(results_diff1, file = paste0(base_dir, "/fc_perm_connectivity1/perm_", perm_id, "/results_diff_scan_", scan_id, ".rds"))
      
    }
  }
}

############
# Load output and save to MATLAB-readable format
############
library(R.matlab)

results_array_task <- array(NA, dim = c(num_edges, num_perm, num_days))
results_array_piecewise <- array(NA, dim = c(num_edges, num_perm, num_days))
results_array_connectivity <- array(NA, dim = c(num_edges, num_perm, num_days))

# Movie vs. Rest
for (perm_id in 1:num_perm) {
  # Load the .rds files for scan_id 1 and 2
  results_scan1 <- readRDS(paste0(base_dir, "/fc_perm/perm_", perm_id, "/results_diff_scan_1.rds"))
  results_scan2 <- readRDS(paste0(base_dir, "/fc_perm/perm_", perm_id, "/results_diff_scan_2.rds"))
  
  # Assign the loaded values to the results_array
  results_array_task[, perm_id, 1] <- results_scan1
  results_array_task[, perm_id, 2] <- results_scan2
}
writeMat(paste0(data_dir, "/isc_heritability/data/solar/fc_perm/results_array_task.mat"), results_array_task = results_array_task)

# Piecewise vs. Anatomical
for (perm_id in 1:num_perm) {
  # Load the .rds files for scan_id 1 and 2
  results_scan1 <- readRDS(paste0(base_dir, "/fc_perm_piecewise1/perm_", perm_id, "/results_diff_scan_1.rds"))
  results_scan2 <- readRDS(paste0(base_dir, "/fc_perm_piecewise1/perm_", perm_id, "/results_diff_scan_2.rds"))
  
  # Assign the loaded values to the results_array
  results_array_piecewise[, perm_id, 1] <- results_scan1
  results_array_piecewise[, perm_id, 2] <- results_scan2
}
writeMat(paste0(data_dir, "/isc_heritability/data/solar/fc_perm/results_array_piecewise.mat"), results_array_piecewise = results_array_piecewise)

# Connectivity vs. Anatomical
for (perm_id in 1:num_perm) {
  # Load the .rds files for scan_id 1 and 2
  results_scan1 <- readRDS(paste0(base_dir, "/fc_perm_connectivity1/perm_", perm_id, "/results_diff_scan_1.rds"))
  results_scan2 <- readRDS(paste0(base_dir, "/fc_perm_connectivity1/perm_", perm_id, "/results_diff_scan_2.rds"))
  
  # Assign the loaded values to the results_array
  results_array_connectivity[, perm_id, 1] <- results_scan1
  results_array_connectivity[, perm_id, 2] <- results_scan2
}

writeMat(paste0(data_dir, "/isc_heritability/data/solar/fc_perm/results_array_connectivity.mat"), results_array_connectivity = results_array_connectivity)