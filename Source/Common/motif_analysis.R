
#' @param positions: List of positions from the positive strand
#' @param input_fasta_obj: FASTA object created by using the 'read.fasta' function from the 'seqinr' library
pos_strand_term_seq_substr <- function(positions, input_fasta_obj) {
  positive_motif <- purrr::map2_chr(positions - 15,
                                      positions,
                                      ~ str_replace_all(toupper((
                                        paste0(c(input_fasta_obj[.x:.y]), collapse = "")
                                      )), "T", "U"))

  return(positive_motif)

}


#' @param positions: List of positions from the negative strand
#' @param input_fasta_obj: FASTA object created by using the 'read.fasta' function from the 'seqinr' library
neg_strand_term_seq_substr <- function(positions, input_fasta_obj) {
  negative_motif <- purrr::map2_chr(positions,
                                    positions + 15,
                                    ~ str_replace_all(as.character(reverseComplement(
                                      DNAString(paste0(c(input_fasta_obj[.x:.y]), collapse = ""))
                                    )), "T", "U"))



  return(negative_motif)
}



#' @param master_table: A table with the strand and position columns.
#' @param input_fasta_obj: FASTA object created by using the 'read.fasta' function from the 'seqinr' library
#' @return A table with the position, strand and sequence that is cut out based the positions
term_seq_substr <-function(master_table, input_fasta_obj) {

  positive_positions <- master_table %>%
    dplyr::filter( strand == "+") %>%
    pull(position)

  negative_positions <-  master_table %>%
    dplyr::filter( strand == "-") %>%
    pull(position)

  positive_motif <- pos_strand_term_seq_substr(positive_positions, input_fasta_obj )

  negative_motif <- neg_strand_term_seq_substr(negative_positions, input_fasta_obj)

  positive_table <- data.frame(positions=positive_positions, sequence=positive_motif, strand="+")
  negative_table <-data.frame(positions=negative_positions, sequence=negative_motif, strand="-")

  full_table <- rbind(positive_table, negative_table)

  return(full_table)
}
