## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
base_dir <- here::here()

data_dir <- file.path( base_dir, "Data" )

graph_dir <- file.path(data_dir, "Bedtools_Count")

#results_dir <- file.path(base_dir, "Results")
results_dir <- "/srv/scratch/treelab/SW_temporary_TermSeq_Results_of_run_pick_peak"
#results_dir <- "/Data"

source_dir <- file.path(base_dir, "Source")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(!require(pacman)) {
  install.packages("pacman")
  libary(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(here)
p_load(UpSetR)

source( file.path( source_dir, "Common", "common_functions.R") )
source( file.path( source_dir, "Common", "group_positions.R") )



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

use_log <- TRUE
get_final_results <- FALSE


compute_platform <- "home"
window_length <- 20   # 2
std_dev_cutoff <- 3 # 7
replicate_cpm_cutoff <- 2
distance_window <- 5
num_conditions_threshold <- 4

# Threshold number of replicates with reads for each genomic position for it to be included in further analysis
replicate_thresholds <- vroom::vroom( file.path( data_dir,
                         "Replicate_Thresholds",
                         "replicate_threshold_b_subtilis.tab" ))

# replicate_thresholds <- vroom::vroom( file.path( data_dir, 
#                          "Replicate_Thresholds", 
#                          "replicate_threshold_s_aureus_jkd6009.tab" ))

column_names_array <- replicate_thresholds %>% distinct(condition) %>% 
  arrange(condition) %>% pull(condition)  

options <- commandArgs(trailingOnly = TRUE)

if(length(options) > 0 ) {
  
  compute_platform <- options[1]
  window_length <- as.numeric(options[2])
  std_dev_cutoff <- as.numeric(options[3])
  
} 

## Parameters for Looping
# cpm_cutoff_list <- c(1, 2, 4, 8)
# distance_window_list <- seq( 5, 20, 5)
cpm_cutoff_list <- c(2, 4)
distance_window_list <- c(5)


if( get_final_results == TRUE) {
  window_length <- 8
  std_dev_cutoff <-  7
  cpm_cutoff_list <- c(1)
  distance_window_list <- c( 5)
  
  if(use_log==TRUE) {
    window_length <- 8
    std_dev_cutoff <-  4
    cpm_cutoff_list <- c(1)
    distance_window_list <- c( 5)    
  }
}

if( use_log == TRUE) {
  
  cpm_cutoff_list <- log2(cpm_cutoff_list)
}

input_parameter_string <-  paste( "win", window_length,  "_",
                                  "sd", std_dev_cutoff,  
                                  sep="" ) 

print(paste( "compute_platform = ", compute_platform, ", ",
             "window_length = ", window_length, ", ",
             "std_dev_cutoff = ", std_dev_cutoff,
             sep=""))



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

upset_counts_dir <- file.path(results_dir, "Sites", "UpSet_Counts")
upset_ratios_dir <- file.path(results_dir, "Sites", "UpSet_Ratios")
selected_params_dir <- file.path(results_dir, "Sites", "Selected_Parameters")


pick_peak_dir <- file.path( results_dir, 
                                  "Sites", 
                                  "Pick_Peak") 

if (use_log == TRUE ) {

  upset_counts_dir <- file.path(results_dir, "Sites", "Log_UpSet_Counts")
  upset_ratios_dir <- file.path(results_dir, "Sites", "Log_UpSet_Ratios")
  selected_params_dir <- file.path(results_dir, "Sites", "Log_Selected_Parameters")

  pick_peak_dir <- file.path( results_dir, 
                                  "Sites", 
                                  "Log_Pick_Peak") 
  
}

create_dir_if_not_exists(upset_counts_dir)
create_dir_if_not_exists(upset_ratios_dir)
create_dir_if_not_exists(selected_params_dir)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cpm_long_tab_file <-  "cpm_long_tbl.tab"

if(use_log == TRUE) {
  
  cpm_long_tab_file <-  "log_cpm_long_tbl.tab"

}

if (! file.exists(file.path( results_dir, 
             "Normalized_Data", 
             "cpm_only_unpivoted.RDS"))) {

cpm_only_unpivoted_temp  <- vroom::vroom(  
  file.path( results_dir, 
             "Normalized_Data", 
             cpm_long_tab_file)) 


cpm_only_unpivoted <- cpm_only_unpivoted_temp %>%
  mutate( Sample = 
            paste( "norm_depth", condition,  replicate,  case_when ( strand == "+" ~ "pos",
                                                                     strand == "-" ~ "neg",
                                                                     TRUE ~ NA_character_), sep="_" )    ) %>%
  arrange(Sample, strand, position)

rm( cpm_only_unpivoted_temp)
gc()

 saveRDS(cpm_only_unpivoted, 
         file.path( results_dir, 
             "Normalized_Data", 
             "cpm_only_unpivoted.RDS"))
} else {
  
  cpm_only_unpivoted <- readRDS( file.path( results_dir, 
             "Normalized_Data", 
             "cpm_only_unpivoted.RDS"))
}




## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## List all the peaks
spikes <- vroom::vroom(file.path( pick_peak_dir, 
                                  paste0("cpm_only_wide_peaks_", input_parameter_string, ".tab")))


## Convert spikes table into a table with position and strand columns
all_positions_temp <- get_all_positions_temp (spikes)

# List of column names 
matrix_column_names_array <-  all_positions_temp %>% distinct(condition) %>% arrange(condition) %>% pull(condition)

## Initialise output table for upset counts data frame
num_expt_cond <- length(matrix_column_names_array)
num_rows_mat <- (2^num_expt_cond)*length(cpm_cutoff_list)*length(distance_window_list)
output_matrix <- as.data.frame( matrix( NA, ncol= num_expt_cond + 3, nrow= num_rows_mat ))
colnames( output_matrix ) <- c( matrix_column_names_array,	"counts", "cpm_cutoff", "distance_window")
rows_count <- 0

## Initialise output table for jaccard ratio data frame
num_rows_ratio <- length(cpm_cutoff_list)*length(distance_window_list)
ratio_matrix <- data.frame( matrix( NA, ncol=3, nrow=num_rows_ratio) )
colnames( ratio_matrix) <- c("cpm_cutoff", "distance_window", "ratio" )
ratio_rows_count <- 0 


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for(replicate_cpm_cutoff in cpm_cutoff_list ) {
  ## Filter by cpm in replicate
  ### This needs loop
  cpm_only_filtered_pos_pivoted <- get_cpm_pivoted_table( cpm_only_unpivoted, 
                                                          replicate_thresholds, 
                                                          replicate_cpm_cutoff) 

  all_positions <- get_all_positions ( all_positions_temp, cpm_only_filtered_pos_pivoted) 
  
  rm(cpm_only_filtered_pos_pivoted)
  gc()
  
  if( get_final_results == TRUE) {
    vroom::vroom_write( all_positions, 
                        file.path( selected_params_dir, paste0("all_positions_and_cpm_", input_parameter_string, 
                                                               "_cpm", replicate_cpm_cutoff, ".tab")) )
  }
  
  if ( nrow(all_positions) == 0 ) {
    next
  }
  
  for (distance_window in distance_window_list) {
    
    print( paste( "cpm_cutoff = ", replicate_cpm_cutoff, 
                  ", distance_window = ", distance_window))
    
    ### This needs loop
    # print("Group all peaks")
    groups_for_all_peaks <- get_group_for_all_peaks( all_positions, distance_window) 
    
    
    if ( nrow( groups_for_all_peaks) == 0 ) {
      next
    }
    
    full_parameter_string <-  paste( "win", window_length,  "_",
                                  "sd", std_dev_cutoff,  "_",
                                  "cpm", replicate_cpm_cutoff, "_",
                                  "gw", distance_window,
                                  sep="" ) 
  
    if( get_final_results == TRUE) {
      vroom::vroom_write( groups_for_all_peaks, 
                          file.path( selected_params_dir, paste0( "groups_for_all_peaks_", 
                                     full_parameter_string, ".tab")))
    }
    
    ## Conditions versus positions
    print( "Conditions versus positions")
    condition_vs_position_clusters <- get_cond_vs_pos_clusters( all_positions, 
                                                                groups_for_all_peaks, 
                                                                replicate_thresholds)
    
    if ( nrow(condition_vs_position_clusters) == 0 ) {
      next
    }
    
    if( get_final_results == TRUE) {
      vroom::vroom_write( condition_vs_position_clusters, 
                          file.path( selected_params_dir, paste0("condition_vs_position_clusters_", 
                                     full_parameter_string, ".tab")) )
    }
    
    
    rm(groups_for_all_peaks)
    gc()
    
    if( get_final_results != TRUE) {
      ## Venn diagram counts 
      print("Venn diagram counts")
      upset_counts <- get_upset_count(condition_vs_position_clusters, column_names_array)
      
      if ( nrow(upset_counts) == 0 ) {
          next
      }
      
      rm(condition_vs_position_clusters)
      gc()
      
      upset_counts_updated <- upset_counts %>%
        mutate( cpm_cutoff = replicate_cpm_cutoff,
                distance_window = distance_window)
      
      output_matrix[(rows_count+1):(rows_count+nrow(upset_counts_updated)), 
                    colnames( upset_counts_updated) ] <- upset_counts_updated
      
      rows_count <- rows_count + nrow(upset_counts_updated)
      
      ## Intersection counts over total counts ratio 
      print( "Get ratio")
      output_ratio <- get_intersection_ratio( upset_counts, column_names_array,
                                              num_conditions_threshold)  
      
      
      if ( length(output_ratio) == 0 ) {
          next
      }
      
      
      ratio_rows_count <- ratio_rows_count + 1
      ratio_matrix[ratio_rows_count, ] <- c(replicate_cpm_cutoff, distance_window, output_ratio )
      
    }
  }
  
    rm(all_positions)
    gc()
}


if( get_final_results != TRUE) {
  
  output_matrix <- output_matrix[1:rows_count,]
  
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if( get_final_results != TRUE) {
  
  vroom::vroom_write( output_matrix, 
                      file.path( upset_counts_dir, 
                                 paste( "upset_counts_", 
                                        input_parameter_string, 
                                        ".tab", sep="")  ) )
  
  vroom::vroom_write( ratio_matrix,
                      file.path( upset_ratios_dir, 
                                 paste( "upset_ratios_", 
                                        input_parameter_string, 
                                        ".tab", sep="") ) )

}


