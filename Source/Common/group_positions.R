
## Group positions that are within a set distance of each other
group_positions <- function( input_table,
                             threshold_distance = 5) {

  separation_dist <- c(threshold_distance+1, input_table$position[-1]  - input_table$position[-length(input_table$position)])

  above_dist_cutoffs <- sapply( separation_dist, function(x){ ifelse(x > threshold_distance, 1, 0 )   } )

  my_groups <- data.frame( group = cumsum( above_dist_cutoffs))

  return( my_groups)

}

# Get the group table min and max read count per position across different experimental conditions
grouped_table <- function( input_table ) {

  all_positions_clustered_summary <- input_table  %>%
    mutate(strand_position =  paste0( strand, position) )  %>%
    group_by ( group ) %>%
    summarise( num_position = n(),
               min_position =  min(position),
               max_position = max(position),
               distance = max(position ) - min(position) + 1,
               position_string = paste0(strand_position, collapse = ", "))  %>%
    ungroup()


  return( all_positions_clustered_summary )
}


## Filter by cpm in replicate
get_cpm_pivoted_table <-  function( cpm_only_unpivoted, replicate_thresholds, replicate_cpm_cutoff) {

  cpm_only_filtered_pos <- cpm_only_unpivoted %>%
    dplyr::filter( norm_depth >  replicate_cpm_cutoff ) %>%
    group_by(  condition, strand, position  ) %>%
    summarise( counts = n()) %>%
    ungroup() %>%
    left_join( replicate_thresholds %>%
               dplyr::select( condition, min_count_num_replicates_threshold),
               by=c("condition")) %>%
    dplyr::filter ( counts >= min_count_num_replicates_threshold )


  cpm_only_filtered_pos_pivoted   <- cpm_only_unpivoted %>%
    inner_join( cpm_only_filtered_pos %>%
                  dplyr::select( -counts),
                by=c( "condition", "strand", "position"))

  return(cpm_only_filtered_pos_pivoted)

}



# Pivot longer the positions
get_all_positions_temp <- function(spikes) {

  all_positions_temp <- spikes %>%
    mutate( position = row_number() ) %>%
    pivot_longer( cols=contains("norm"),
                  names_to="Sample",
                  values_to="Is_Spike") %>%
    dplyr::filter( Is_Spike == TRUE) %>%
    mutate( condition = str_replace( Sample, "norm_depth_(.*?)_.*", "\\1")  ) %>%
    mutate( replicate = str_replace( Sample, "norm_depth_.*?_(\\d+)_(pos|neg)", "\\1")  ) %>%
    mutate( strand = str_replace( Sample, "norm_depth_.*?_\\d+_(pos|neg)", "\\1")  )  %>%
    mutate( strand = case_when (  strand == "neg"~ "-",
                                  strand == "pos"~ "+",
                                  TRUE ~ NA_character_))

  return( all_positions_temp)

}



# Join positions with their CPM data
get_all_positions <- function( all_positions_temp, cpm_only_filtered_pos_pivoted) {
  all_positions <- all_positions_temp %>%
    mutate( replicate = as.double(replicate)) %>%
    inner_join( cpm_only_filtered_pos_pivoted ,
                by=c("position",
                     "Sample",
                     "condition",
                     "replicate",
                     "strand") )  %>%
    arrange(Sample, strand, position)

  return( all_positions)

}


# Group the peaks that are within a distance window of each other
get_group_for_all_peaks <- function( all_positions, distance_window) {

  groups_for_all_peaks_temp <-   all_positions %>%
    distinct( strand, position) %>%
    arrange( strand, position) %>%
    group_by(strand ) %>%
    nest() %>%
    mutate( data2 = purrr::map( data, ~group_positions(., threshold_distance = distance_window))) %>%
    unnest( c(data, data2 ) ) %>%
    ungroup()

  max_neg_group_number <- groups_for_all_peaks_temp %>%
    dplyr::filter( strand == "-" ) %>%
    summarise( max( group)) %>%
    pull()

  # Merge group number for positive and negative strand
  groups_for_all_peaks <- groups_for_all_peaks_temp %>%
    dplyr::mutate(  group = case_when( strand == "-" ~ group,
                                       strand == "+" ~ group + max_neg_group_number,
                                       TRUE ~ NA_real_ ))

  return( groups_for_all_peaks )
}


# Map groups of peaks to different experimental conditions
# Each peak must be in the same group for two or more biological replicates
get_cond_vs_pos_clusters <- function(all_positions, groups_for_all_peaks, replicate_thresholds) {

  condition_vs_position_clusters <- all_positions %>%
    distinct( condition, position, replicate, strand) %>%
    left_join(groups_for_all_peaks, by=c("position", "strand") ) %>%
    group_by(condition,  group) %>%
    summarise( num_replicates = n()) %>%
    ungroup()  %>%
    left_join( replicate_thresholds %>%
                dplyr::select( condition, peaks_num_replicates_threshold),
               by=c("condition")) %>%
    dplyr::filter ( num_replicates >= peaks_num_replicates_threshold )%>%
    dplyr::select(-num_replicates) %>%
    distinct() %>%
    mutate( is_present = 1) %>%
    arrange( group, condition ) %>%
    pivot_wider( id_cols=c(group),
                 names_from=c(condition),
                 values_from = c(is_present),
                 values_fill = list( is_present  = 0   )  )

  return( condition_vs_position_clusters )
}


## Use across() function
## how many peaks were found amoung this list of conditions
get_upset_count <- function(condition_vs_position_clusters, column_names_array) {


  shared_col_names <- intersect(colnames(condition_vs_position_clusters ), column_names_array)

  ## Upset counts
  upset_counts  <- condition_vs_position_clusters %>%
    group_by( across( one_of( shared_col_names)) ) %>%
    summarise( counts = n()) %>%
    ungroup()

  return(upset_counts )
}

## Get jaccard ratio. Num peaks found in all conditions versus total number of peaks
get_intersection_ratio <- function( upset_counts, column_names_array, num_conditions_threshold) {

#   intersection_counts <- upset_counts %>%
#     dplyr::filter( cont  == 1 &
#                      linz  == 1 &
#                      vanco == 1 &
#                      tige == 1 ) %>%
#     summarise( my_total_counts =  sum(counts)) %>%
#     dplyr::pull(my_total_counts)

  shared_col_names <- intersect(colnames(upset_counts ), column_names_array)

  get_is_present <- upset_counts %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(cols= all_of(shared_col_names),
                 names_to = "treatment",
                 values_to = "is_present") %>%
    dplyr::filter( is_present == 1)

  common_to_all <- get_is_present %>%
    group_by(row_id) %>%
    summarise( num_expt = n()) %>%
    ungroup() %>%
    dplyr::filter( num_expt >= num_conditions_threshold ) %>%
    dplyr::select(-num_expt)

  intersection_counts <- get_is_present %>%
    dplyr::inner_join( common_to_all, by=c("row_id")) %>%
    distinct( counts, row_id) %>%
    pull(counts)


  total_counts <- upset_counts %>%
    summarise( my_total_counts =  sum(counts)) %>%
    pull(my_total_counts)

  return( intersection_counts/total_counts )


}
