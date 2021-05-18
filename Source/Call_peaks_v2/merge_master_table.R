## -------------------------------------------------------------------------------------------------------------------------------------------------------------
base_dir <- here::here()

data_dir <- file.path( base_dir, "Data" )

graph_dir <- file.path(data_dir, "Bedtools_Count")

# results_dir <- file.path(base_dir, "Results")
results_dir <- "/srv/scratch/treelab/SW_temporary_TermSeq_Results_of_run_pick_peak"
# results_dir <- file.path(base_dir, "Results", "Results_B_subtilis")

source_dir <- file.path(base_dir, "Source")



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

if(!require(pacman)) {
  install.packages("pacman")
  libary(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(here)
p_load(UpSetR)
p_load(magrittr)

source( file.path( source_dir, "Common", "common_functions.R") )
source( file.path( source_dir, "Common", "group_positions.R") )



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
  use_log <- TRUE
  window_length <- 8
  std_dev_cutoff <-  7
  replicate_cpm_cutoff <- 1
  distance_window <- 5
  
  
if(use_log==TRUE) {
  
  window_length <- 8
  std_dev_cutoff <-  3
  replicate_cpm_cutoff <- 1
  distance_window <- 5
  
  replicate_cpm_cutoff <- log2(replicate_cpm_cutoff)
}  

input_parameter_string <-  paste( "win", window_length,  "_",
                                  "sd", std_dev_cutoff,  
                                  sep="" ) 

full_parameter_string  <-  paste( "win", window_length,  "_",
                                  "sd", std_dev_cutoff,  "_",
                                  "cpm", replicate_cpm_cutoff, "_",
                                  "gw", distance_window,
                                  sep="" ) 

print(paste(  "window_length = ", window_length, ", ",
             "std_dev_cutoff = ", std_dev_cutoff,
             sep=""))



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

selected_params_dir <- file.path(results_dir, "Sites", "Selected_Parameters")

if (use_log == TRUE ) {
  selected_params_dir <- file.path(results_dir, "Sites", "Log_Selected_Parameters")
}

create_dir_if_not_exists(selected_params_dir)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
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





## -------------------------------------------------------------------------------------------------------------------------------------------------------------

all_positions <- vroom::vroom(
  file.path( selected_params_dir, paste0("all_positions_and_cpm_", input_parameter_string, 
                                         "_cpm", replicate_cpm_cutoff, ".tab")) )

groups_for_all_peaks <- vroom::vroom( 
  file.path( selected_params_dir, paste0( "groups_for_all_peaks_", 
                                          full_parameter_string, ".tab")))

condition_vs_position_clusters <- vroom::vroom( 
  file.path( selected_params_dir, paste0("condition_vs_position_clusters_", 
                                     full_parameter_string, ".tab")) )



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
all_positions_cln <- all_positions %>%
  dplyr::select( -Sample, -Is_Spike) %>%
  mutate( sample = paste( condition, replicate, sep="_") ) %>%
  dplyr::select( -condition, -replicate) %>%
  pivot_wider( id_cols=c(position, strand),
               names_from="sample",
               values_from=c("depth", "norm_depth"))


rm( all_positions)
gc()



## -------------------------------------------------------------------------------------------------------------------------------------------------------------


log_cpm_wider_tbl_temp <-  cpm_only_unpivoted %>%
  inner_join( all_positions_cln %>% 
                dplyr::select(position, strand),
              by=c("position", "strand") ) 


rm(cpm_only_unpivoted)
gc()


cpm_pos_max <- log_cpm_wider_tbl_temp %>%
  group_by( position, strand ) %>%
  summarise( max_norm_depth = max(norm_depth)) %>%
  ungroup()
  


## -------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Find position with maximum posible read depth per group of peaks ")

max_depth_per_group <-  groups_for_all_peaks %>%
  inner_join( all_positions_cln %>% 
                distinct(position, strand), by=c("position", "strand")) %>%
  left_join( cpm_pos_max, by=c( "position", "strand") ) %>%
  group_by( group,  strand) %>%
  summarise( max_norm_depth = max(max_norm_depth)) %>%
  ungroup() %>%
  mutate( is_max_position = 1)

## If there are ties in the max normalized depth, choose the earliest position (smallest position number)
max_position_per_group_temp <-  groups_for_all_peaks %>%
  left_join( cpm_pos_max, by=c( "position", "strand") ) %>%
  left_join (max_depth_per_group, by = c("group",  "strand", "max_norm_depth")) %>%
  dplyr::filter( is_max_position == 1) 

  

max_position_per_group <- max_position_per_group_temp %>%
  group_by(strand,  group, max_norm_depth, is_max_position ) %>%
  summarise( position = min(position)) %>%
  ungroup()




## -------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("log cpm wider step")


log_cpm_wider_tbl <- log_cpm_wider_tbl_temp %>%
  dplyr::select(-Sample  ) %>%
  dplyr::rename( orig_depth = "depth",
                 orig_norm_depth = "norm_depth") %>%
  pivot_wider( id_cols=c( position, strand ),
               names_from=c( condition, replicate),
               values_from = c(orig_depth, orig_norm_depth))


rm(log_cpm_wider_tbl_temp)
gc()



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
cat("Merge master table")

master_table_before_filter <- all_positions_cln %>%
  left_join( groups_for_all_peaks, by=c('position', 'strand') ) %>%
  left_join( condition_vs_position_clusters, by="group") %>%
  left_join( max_position_per_group, by=c("strand", "position", "group")) %>%
  left_join( log_cpm_wider_tbl, by=c("strand", "position")) %>%
  mutate(  is_max_position = ifelse( is.na(is_max_position), 0, 1  )) %>%
  dplyr::rename( prelim_cont = "cont",
                 prelim_linz = "linz",
                 prelim_tige = "tige",
                 prelim_vanco = "vanco") %>%
  dplyr::mutate( cont = case_when(  orig_norm_depth_cont_1 <= 0 &  orig_norm_depth_cont_2 <= 0 ~ 0,
                                   prelim_cont == 1 ~ 1,
                                   TRUE ~ 0 ),
                 linz = case_when( orig_norm_depth_linz_1 <= 0 & orig_norm_depth_linz_2  <= 0 ~ 0,
                                   prelim_linz == 1 ~ 1,
                                   TRUE ~ 0  ),
                 tige = case_when( orig_norm_depth_tige_1 <= 0 & orig_norm_depth_tige_2  <= 0 ~ 0,
                                   prelim_tige == 1 ~ 1,
                                   TRUE ~ 0  ),
                 vanco = case_when( orig_norm_depth_vanco_1 <= 0 &  orig_norm_depth_vanco_2 <= 0 ~ 0,
                                    prelim_vanco == 1 ~ 1,
                                   TRUE ~ 0  )  ) %>%
  dplyr::select( group, position, strand, is_max_position, cont, linz, tige, vanco, contains("depth"),
                 prelim_cont, 
                 prelim_linz,
                 prelim_tige,
                 prelim_vanco ) %>%
  arrange(group, position, strand)

master_table <-  master_table_before_filter %>%
  dplyr::filter( is_max_position == 1)  %>%
  dplyr::select(-is_max_position)



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

vroom::vroom_write( master_table_before_filter, 
  file.path( selected_params_dir, 
             paste0("master_table_unfiltered_", 
                    full_parameter_string, ".tab")))



## -------------------------------------------------------------------------------------------------------------------------------------------------------------

vroom::vroom_write( master_table, 
  file.path( selected_params_dir, 
             paste0("master_table_", 
                    full_parameter_string, ".tab")))


## -------------------------------------------------------------------------------------------------------------------------------------------------------------

master_table_upset <- master_table %>%
  mutate( position_strand = paste(position, strand, sep="_")) %>%
  dplyr::select( position_strand, cont, linz, tige, vanco) %>%
  column_to_rownames("position_strand") 

upset(master_table_upset, order.by = "freq")

tiff( file.path(selected_params_dir, "upset_freq_master_table.tiff" ))
upset(master_table_upset, order.by = "freq")
dev.off()

upset(master_table_upset, order.by = "degree")

pdf( file.path(selected_params_dir, "upset_degree_master_table.pdf" ))
upset(master_table_upset, order.by = "degree")
dev.off()





## -------------------------------------------------------------------------------------------------------------------------------------------------------------

cont_term_pos <- master_table %>%
  dplyr::filter( strand == "+" & cont == 1) %>%
  dplyr::select( position, norm_depth_cont_1, norm_depth_cont_2 ) %>%
  dplyr::mutate( max_norm_depth = case_when(norm_depth_cont_1 > norm_depth_cont_2 ~ norm_depth_cont_1, 
                                            TRUE ~ norm_depth_cont_2 )) %>%
  dplyr::select( -norm_depth_cont_1, -norm_depth_cont_2 ) %>%
  #dplyr::filter( !is.na(max_norm_depth) ) %>%
  set_colnames(  c("variableStep", "chrom=Saa_6009_C1") )


cont_term_neg <- master_table %>%
  dplyr::filter( strand == "-" & cont == 1) %>%
  dplyr::select( position, norm_depth_cont_1, norm_depth_cont_2 ) %>%
  dplyr::mutate( max_norm_depth = case_when(norm_depth_cont_1 > norm_depth_cont_2 ~ norm_depth_cont_1, 
                                            TRUE ~ norm_depth_cont_2 )) %>%
  dplyr::select( -norm_depth_cont_1, -norm_depth_cont_2 ) %>%
  #dplyr::filter( !is.na(max_norm_depth) ) %>%
  set_colnames(  c("variableStep", "chrom=Saa_6009_C1") )


nrow(cont_term_pos) + nrow(cont_term_neg)

vroom::vroom_write( cont_term_pos, file.path(selected_params_dir, "Wiggle_Format", "cont_term_pos.wig"  ))

vroom::vroom_write( cont_term_neg, file.path(selected_params_dir, "Wiggle_Format", "cont_term_neg.wig"  ))



## -------------------------------------------------------------------------------------------------------------------------------------------------------------
all_term_pos <- master_table %>%
  dplyr::filter( strand == "+") %>%
  dplyr::select( position, max_norm_depth ) %>%
  set_colnames(  c("variableStep", "chrom=Saa_6009_C1") )


all_term_neg <- master_table %>%
  dplyr::filter( strand == "-") %>%
  dplyr::select( position, max_norm_depth ) %>%
  set_colnames(  c("variableStep", "chrom=Saa_6009_C1") )


vroom::vroom_write( all_term_pos, file.path(selected_params_dir, "Wiggle_Format", "all_term_pos.wig"  ))

vroom::vroom_write( all_term_neg, file.path(selected_params_dir, "Wiggle_Format", "all_term_neg.wig"  ))

# /home/ignatius/PostDoc/2020/TermSeq2020/Results/Sites/Log_Selected_Parameters/Wiggle_Format

