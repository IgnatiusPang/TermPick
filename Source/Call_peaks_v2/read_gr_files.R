## ------------------------------------------------------------------------

if(!require(pacman)) {
  install.packages("pacman")
  libary(pacman)
}

p_load(tidyverse)
p_load(vroom)
p_load(here)
p_load(devtools)

if(!require(multidplyr)) {
  devtools::install_github("tidyverse/multidplyr")
  library(multidplyr)
}





## ------------------------------------------------------------------------

create_dir_if_not_exists <- function(file_path, mode = "0777") { 
  
  ### Create directory recursively if it doesn't exist
  if (! file.exists(file_path)){
    
    dir.create(file_path, showWarnings = TRUE, recursive = TRUE, mode = mode)
    
  }
  
}



## ------------------------------------------------------------------------
base_dir <- here::here()

data_dir <- file.path( base_dir, "Data" )

graph_dir <- file.path(data_dir, "Bedtools_Count")

results_dir <- file.path(base_dir, "Results")

create_dir_if_not_exists(file.path(results_dir, "Normalized_Data"))

source_dir <- 




## ------------------------------------------------------------------------
cluster <- new_cluster(6)

read_count_threshold <- 4
read_cpm_threshold <- 4



## ------------------------------------------------------------------------
list_of_files <-  Sys.glob( file.path( graph_dir, "*.gr") ) 

list_of_files


## ------------------------------------------------------------------------

term_pos_tbl  <- vroom::vroom( list_of_files, 
                               id="filename", 
                               col_names=c("genome_id", "position", "depth") )




## ------------------------------------------------------------------------
term_pos_cleaned <- term_pos_tbl %>%
      mutate(  cleaned_f_name = str_replace ( filename, paste0( graph_dir, "/" ), ""    ))  %>%
      dplyr::select( cleaned_f_name, position, depth)



## ------------------------------------------------------------------------
term_pos_condition <- term_pos_cleaned %>%
  mutate( strand = case_when( str_detect(cleaned_f_name, "minus")  ~ "-",
                              str_detect( cleaned_f_name, "plus" ) ~ "+",
          TRUE ~ NA_character_ ) ) %>%
  mutate( condition = str_replace( cleaned_f_name, "Staph_termseq_(.*?)(\\d?)_(.*)", "\\1"   ),
          replicate = str_replace( cleaned_f_name, "Staph_termseq_(.*?)(\\d?)_(.*)", "\\2"   )
          ) %>%
  # Clean replicate labels, in which replicate 3 and 4 are replaced with replicate 1 and 2 accordingly
  mutate( replicate  = case_when(replicate == "3" ~ "1",
                                 replicate == "4" ~ "2",
                                 TRUE ~ replicate )) %>%
  dplyr::select( condition, replicate, position, strand, depth)  %>%
  arrange(condition, replicate, position, strand )

# Staph_termseq_vanco1_minus_strand_single_3end_base_only.gr





## ------------------------------------------------------------------------

vroom::vroom_write(term_pos_condition, path=file.path(results_dir, "Staph_termseq_3end_base_only.tab"),
                   delim = "\t")



## ------------------------------------------------------------------------

total_depth_per_sample <- term_pos_condition %>%
  group_by(condition, replicate) %>%
  summarise( total_depth = sum(depth)) %>%
  ungroup()

total_depth_per_sample


## ------------------------------------------------------------------------
cpm_long_tbl <- term_pos_condition %>%
  left_join( total_depth_per_sample, by =c( "condition", "replicate")) %>%
  # Calculate the counts per million value 
  mutate( norm_depth =  ( (depth + 1) / total_depth * 10^6 )  ) %>%
  dplyr::select(-total_depth) %>%
  as_tibble()

cpm_long_tbl


plot(density(cpm_long_tbl$norm_depth))


plot(density(2^cpm_long_tbl$norm_depth))




## ------------------------------------------------------------------------
vroom::vroom_write( cpm_long_tbl, 
                    file.path( results_dir, 
                               "Normalized_Data", 
                               "cpm_long_tbl.tab"))


## ------------------------------------------------------------------------
cpm_wider_tbl <-  cpm_long_tbl  %>%
  dplyr::mutate(  strand = case_when ( strand == "+" ~ "pos",
                                       strand == "-" ~ "neg",
                                       TRUE ~ NA_character_)) %>%
  pivot_wider( id_cols=c( position ),
               names_from=c( condition, replicate, strand),
               values_from = c(depth, norm_depth))




## ------------------------------------------------------------------------
vroom::vroom_write( cpm_wider_tbl, 
                    file.path( results_dir, 
                               "Normalized_Data", 
                               "cpm_wider_tbl.tab"))


## ------------------------------------------------------------------------

log_cpm_long_tbl <- cpm_long_tbl %>%
  mutate( norm_depth = log2(norm_depth) ) 


vroom::vroom_write( log_cpm_long_tbl, 
                    file.path( results_dir, 
                               "Normalized_Data", 
                               "log_cpm_long_tbl.tab"))



## ------------------------------------------------------------------------
log_cpm_wider_tbl <-  log_cpm_long_tbl %>%
  dplyr::mutate(  strand = case_when ( strand == "+" ~ "pos",
                                       strand == "-" ~ "neg",
                                       TRUE ~ NA_character_)) %>%
  pivot_wider( id_cols=c( position ),
               names_from=c( condition, replicate, strand),
               values_from = c(depth, norm_depth))




## ------------------------------------------------------------------------
vroom::vroom_write( log_cpm_wider_tbl, 
                    file.path( results_dir, 
                               "Normalized_Data", 
                               "log_cpm_wider_tbl.tab"))

